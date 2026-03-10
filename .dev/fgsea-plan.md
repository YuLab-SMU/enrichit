# fgsea 算法高保真复现计划（enrichit）

## 1. 目标与验收标准

### 1.1 目标
- 在 `enrichit` 中实现与 `fgsea` 一致的预排序 GSEA 多层自适应蒙特卡洛算法。
- 支持 `scoreType = "std" | "pos" | "neg"`，并在极小 p 值场景保持数值稳定。
- 提供可复现实验与系统评测，优先保证统计准确性，再优化性能。

### 1.2 验收标准（必须全部满足）
- 统计一致性：相同输入与随机种子下，`enrichit` 与 `fgsea` 的 `ES/NES/pval/log2err` 高度一致。
- 尾部精度：极小 p 值区域（例如 `<1e-12`）保持可解释误差，且与 `fgsea` 同数量级。
- 校准性：null 仿真下 p 值分布接近均匀，BH-FDR 控制无明显失真。
- 回归稳定：新增测试覆盖核心路径，R 包检查通过。

## 2. fgsea 源码阅读结论（复现基线）

基于 `fgsea` 源码（`R/fgseaMultilevel.R`, `R/fgsea.R`, `src/fgseaMultilevel.cpp`）可归纳出关键机制：

1. 两阶段策略  
   - 阶段 A：`fgseaSimpleImpl` 用 `nPermSimple` 粗估计。  
   - 阶段 B：仅对 `multError < simpleError` 的 pathways 进入 multilevel 精算。

2. `std` 模式概率归一化  
   - `denomProb = (modeFraction + 1) / (nPermSimple + 1)`  
   - 最终 `pval = min(1, cppMPval / denomProb)`。

3. 误差模型  
   - `simpleError` 使用 Beta 区间（`qbeta`）与 crude estimator。  
   - `multilevelError(pval, sampleSize)` 用于 multilevel 结果的 `log2err`。

4. 健壮性处理  
   - `modeFraction < 10` 时置 NA。  
   - `pval < eps` 截断为 `eps`。  
   - `sampleSize` 强制奇数且至少为 3。

## 3. 当前 enrichit 的关键差距

1. simple 与 multilevel 的切换逻辑未完全按 `fgsea` 复刻。  
2. stats 整数缩放策略与 `fgsea` 不一致，可能影响尾部精度。  
3. 边界行为（`modeFraction < 10`、`eps`、一侧检验）仍需逐项对齐。  
4. 现有测试偏“可运行”，缺少“高精度一致性与校准性”验证。

## 4. 实施规划（具体到代码层）

### Phase 0：建立对照基线
- 锁定 `fgsea` 版本，固定对照环境与随机种子。
- 建立统一对照脚本，输出 `ES/NES/pval/log2err` 差异表。

### Phase 1：算法行为强一致改造
- 在 `R/gsea.R` 增加与 `fgsea` 对齐的预处理：有限值检查、排序、`abs(stats)^exponent`、整数缩放。
- 改造 `src/gsea_multilevel.cpp`：接收已缩放整数 ranks，去除不一致的内部缩放。
- 完整复刻 `simpleError` 与 `multError` 比较逻辑，只将必要通路送入 multilevel。
- 复刻 `std` 模式 `denomProb` 归一化、`isCpGeHalf` 处理、`eps` 截断与 `log2err` 规则。

### Phase 2：边界行为与 API 对齐
- 对齐 `scoreType = std/pos/neg` 的一侧与双侧语义。
- 对齐 `modeFraction < 10` 的 NA/warning 行为。
- 保证输出列语义与 `gseaResult` 转换稳定。

### Phase 3：评测体系（准确性优先）
- 新增高精度 parity 测试、尾部精度测试、null 校准测试。
- 形成可复现报告与验收阈值判断。

## 5. 详细实现步骤（可执行顺序）

1. 新建对照脚本 `analysis/compare_with_fgsea.R`，在同输入同种子下跑双实现。  
2. 重构 `gsea` 预处理流程，将 stats 变换与缩放统一到 R 侧。  
3. 调整 `gsea_multilevel_cpp` 接口与内部流程，移除“待确认”逻辑分支。  
4. 按 `fgsea` 完整实现 `simpleError vs multError` 路由。  
5. 对齐 `denomProb`、`cppMPval`、`isCpGeHalf`、`log2err`、`eps` 行为。  
6. 补齐边界条件：全正/全负 stats、ties、极小 pathway、极大 pathway、空交集。  
7. 新增测试并设定硬阈值。  
8. 运行完整检查并固化评测报告。

## 6. 评测设计（围绕准确性）

### 6.1 数值一致性评测
- 数据：`fgsea` 自带 `examplePathways/exampleRanks` + 自定义合成数据。
- 指标：
  - `ES` 绝对误差
  - `NES` 相对误差
  - `-log10(pval)` 相关性
  - top-k 通路重叠率

建议阈值：
- `cor(-log10(pval_enrichit), -log10(pval_fgsea)) >= 0.995`
- `NES` 符号一致率 100%
- top-50 overlap >= 95%

### 6.2 尾部精度评测
- 构造强富集通路，重点检查 `<1e-8`, `<1e-12` 区域。
- 对比 `pval` 数量级与 `log2err` 变化趋势。

### 6.3 校准性评测
- null 仿真（随机 pathways + 随机 ranks，多批次重复）。
- 指标：
  - QQ 图接近对角线
  - KS 检验不过度偏离均匀分布
  - BH-FDR 实测接近标称值

### 6.4 鲁棒性评测
- ties 比例从 0% 到 30%
- pathway size 跨度覆盖 min/max 边界
- `scoreType` 三模式全覆盖

## 7. 里程碑

- M1（1-2 天）：完成对照基线脚本与数据固化。  
- M2（2-4 天）：完成 simple/multilevel 核心逻辑对齐。  
- M3（2-3 天）：完成边界行为与 API 对齐。  
- M4（2-4 天）：完成系统评测并生成报告。  
- M5（1 天）：收敛参数建议并完成质量门禁。  

## 8. 质量门禁（必须执行）

- `Rscript -e 'source("tests/testthat.R")'`
- `Rscript -e 'devtools::check()'`

若门禁失败，优先修复准确性问题，再处理性能优化。

## 9. 最终交付物

- 代码改造：`R/gsea.R`, `src/gsea_multilevel.cpp/.h`（及必要 util）。  
- 测试：  
  - `tests/testthat/test-gsea-multilevel-parity.R`  
  - `tests/testthat/test-gsea-tail.R`  
  - `tests/testthat/test-gsea-calibration.R`  
- 分析报告：  
  - `analysis/reports/fgsea_parity_summary.csv`  
  - `analysis/reports/fgsea_tail_accuracy.csv`  
  - `analysis/reports/fgsea_calibration_summary.csv`  

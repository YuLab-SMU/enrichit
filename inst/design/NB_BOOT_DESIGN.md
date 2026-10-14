# NB-native 集合检验设计（`dsea_nb`）—— PAPER_PLAN §10 攻点 2

> 本文件最初写在被误删的 `enrichit/inst/design/` 下，此处为重建后的版本。
> 状态：**设计已定，实现未开始**。落地前须先读完 §4 的证伪判据并同意其可执行性。

## 0. 动机（合法 vs 非法）

**合法**（消费 DESeq2 建模特有的信息）：
- NB 的**偏度**：`Var = μ + αμ²` 使低表达基因的零膨胀与右偏无法由 log-CPM + 高斯
  残差模型表达；voom/delta 权重只给方差标度，不给三阶矩。
- **收缩的可信度**：`mcols(x)$dispersion` 是 EB 收缩后的 α̂，带 `fitDF`/`priorDF`；
  同一 α̂ 在不同先验强度下证据量不同，正态近似看不见这一层。

**非法**（不要写进论文）：
- ❌ "修复组成偏差"——实验六 B 已证伪：对照基因归一化后偏差仍在，且我们的方法与
  limma 系同样受影响。本检验**不声称**解决组成偏差。

## 1. 统计量

主选（候选 B，fry 形式）：

    t_S = mean_k( Eff[S, k=1] ) / sqrt( mean_{k>=2} mean_k( Eff[S,k]^2 ) )

分子为对比坐标的集内均值，分母为残差坐标平方的均方根。与 `st_fry_cpp` 的
`M[1]/sqrt(mean(M[-1]^2))` 同构，可直接复用 `st_effects_cpp` 的结果。

敏感性：候选 A（camera 式加权均值 / 相关校正方差）作为对照实现，仅在结论一致时采用主选。

## 2. 零分布：NB 参数化 bootstrap

    c*_{g,b} ~ NB( μ̂_{g,b}, α̂_g )        # 保留每样本各自的 μ̂，不加真实效应
    → 走与观测值完全相同的流水线（size factor、E、权重、QR、统计量）
    → p = (1 + #{ |t*_S| >= |t_S| }) / (B + 1)      # 双侧

**三个必须写死的岔路口**：
1. **size factor**：固定为观测 `sizeFactors(x)`，不从重采样的 counts 重算
   （重算会把组成效应二次混入零分布）。
2. **μ̂**：用 DESeq2 拟合的 `assay(x,"mu")`（MLE），不重新拟合模型。
3. **α̂**：用 `mcols(x)$dispersion`（EB 收缩后），不用未收缩的 gene-wise 估计。

`B = 1000` 默认；每个集合独立计数；**必须记录种子**并在 golden 里锁一份。

## 3. 输入契约

    dsea_nb(x, gene_sets, contrast = NULL, B = 1000, seed = NULL, ...)

- `x`：已 `DESeq()` 的 `DESeqDataSet`（需要 `mu` 与 `dispersion`）
- `DseaData` 需新增字段：`counts`（原始计数）、`mu`、`alpha`
- 统计量与 effects 复用 `st_effects_cpp`；新增 `st_nb_boot_cpp` 只做重采样循环
- **与 `dsea()` 并行**：不改 `dsea()` 的任何分支，不动它的 parity。

## 4. 校准实验与证伪判据

**判据（任一命中即放弃该检验，按 PAPER_PLAN 的 no-go 纪律）**：
- 零集 α 的 95% CI 与 `fry` 的 CI **重叠** → 说明没有增量信息，废；
- α > 0.15 → 反保守，废；
- 单次运行成本 > `fry` × 50 → 工程上不可用，废；
- 小集（NGenes < 5）p 值过度离散（`B=1000` 下分辨率 1e-3）→ 至少在文档里明确警告。

**跑什么**：
1. PAPER_PLAN 四难场景（低表达 / 高离散 / 小集 / 不平衡设计）上的 α 校准，对照 `fry`；
2. airway 置换零分布上的 α，对照现有 `fry`/`roast` 的 0.124–0.129；
3. 若 α 明显低于 fry 且贴近名义水平，才值得继续谈"定制建模有优势"。

## 5. 交付物

- `enrichit/R/dsea_nb.R`
- `enrichit/src/set_test.cpp`：新增 `st_nb_boot_cpp`
- `DseaData` 增 `counts` / `mu` / `alpha` 三个字段（`as_dsea()` 填充）
- `dsea-benchmark/scripts/exp_nb_calib.R`（四难场景校准）
- `dsea-benchmark/scripts/exp_nb_perm.R`（airway 置换对照）

## 6. 边界

- **不做** camera / roast 的 NB 版本：camera 的 competitive 假设与置换式 bootstrap
  概念冲突，roast 本身已是置换检验，加 NB 零分布是重复劳动。
- **不做** NB 精确检验：参数化 bootstrap 已足够，精确版成本不可控。
- 入口必须是**并行**的：`dsea()` 行为一字不改，`dsea_nb()` 是新增函数。

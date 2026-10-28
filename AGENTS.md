# enrichit — Agent 工作说明

适用于 `enrichit/` 包；工作区级通用规则见 `~/.dsh/AGENTS.md`。

## 1. 定位

`clusterProfiler` 家族的**统计引擎层**：ORA、GSEA（multilevel/permutation/adaptive）、
`dsea()`、加权富集、NSEA/MNSEA 网络传播、多组学融合。算法在 C++（`Rcpp`/`RcppEigen`），
R 层管输入校验、结果对象组装、GSON 适配。上层是 `clusterProfiler`/`enrichplot`/`DOSE` 等，
**改这里会波及整个家族**。

## 2. 测试边界

原则：**谁定义契约，谁负责测试。**

| 内容 | 在本包测 |
|---|---|
| 统计核心：p 值、ES/NES、leading edge、multilevel 守卫、`dsea()` 的对比度与排序 | ✅ 必测 |
| 结果对象与过滤：`get_enriched()`、`qvalueCutoff`/`pvalueCutoff` 语义、转换器 | ✅ 必测 |
| 输入校验与退化输入：全零统计量、`eps` 越界、覆盖全列表的集合、非有限值 | ✅ 必测 |
| `setReadable()` / `EXTID2NAME()` | ❌ 在 clusterProfiler 测 |

`setReadable()` 的语义由上层决定（`clusterProfiler` 再导出它），且要测它就得在 `Suggests`
里挂 100+ MB 的 `org.Hs.eg.db`；在本包能写出的断言只有"`AnnotationDbi::select()` 被调用了"，
真正有价值的（与 `enrichGO()` 等入口的配合）只有上层写得出来。**除非语义不再由上层决定，
或有了不依赖大注释包的构造方式，否则不要在本包补测。**

## 3. 验证

包根有 `Makefile`，**先按它走**：`make check`（全量）、`make rd`（文档）。

实测代价：`load_all()` **6 s** / 全套 `test_dir()` **55 s** / `make check` **2 min**。
2 分钟不算贵，**收口直接 `make check`**，别照搬 sclet（25 min）那套层层设防。

```bash
export R_LIBS_USER=~/Rlib    # 本机必须，否则依赖被误判缺失
```

**判据 `0 errors / 0 warnings`。** 唯一已知 NOTE 是 `-mno-omit-leaf-frame-pointer`，
来自发行版 R 自己的 `CFLAGS`（包内无 `src/Makevars`），不要去"修"。

改 C++ 后必须重新 `load_all()`（首次编译约 1 分钟），否则跑的还是旧 `.so`。

退化输入上包会发诊断（`eps` 到下限、统计量全正、q-value 估不出等），这些是**预期行为**：
用 `suppressWarnings()`/`suppressMessages()` 或 `expect_warning()` 显式声明，
否则 `make check` 会因 WARN>0 失败，真正意外的告警也被淹没。

### 断言必须做变异检查，否则不算写完

绿不等于有效：断言可能**根本没测到目标分支**，或**被另一个机制的报错顺带满足**。
唯一判据是**把被测代码改坏，测试必须红**（先 `cp` 备份，改完跑单文件，再还原）。

实测过的四类失效写法：

- **`for` 遍历 `c(...)` 做非法值用例**：混合向量被强制成 `character`，每轮只走"类型不对"分支，
  越界用例（如 `eps = 2`）从未执行 → **用 `list(...)`**
- **`expect_error` 的 regexp 太宽**：两个守卫的报错都含同一句话时，删掉前一道测试照样绿
  → regexp 要能唯一确定是**哪一层**拦下的
- **`expect_warning(..., ".")`**：任何无关告警都通过，等于没写 → `capture_warnings()` + 逐条 `expect_match()`
- **依赖随机结果的具体值**：先算出实际值再写死，或改断言**结构性质**（如"p 值必落在
  `1/(n+1)` 网格上"）。注意 `%in%` 优先级高于 `/`，`x %in% (1:4)/4` 是 `(x %in% 1:4)/4`，**必须加括号**

正确姿势：期望值要能**由输入手算**（spec），不要从被测对象的输出或实现内部同一数据源反推（那是自证）；
**输入必须完全指定**（`head(keys(odb,"ENTREZID"),120)` 不是受控输入）；
只有 mock 才能到达的分支用 `with_mocked_bindings(..., .package = "enrichit")`，别用"不可达所以不测"糊过去。

## 4. 变更约定

- **版本号**：开发期「正式号 + `.NNN`」每次 commit +1；正式号由作者定。
- **NEWS.md**：顶部小节即当前版本；新条目置顶写在该节内；**已发布小节一个字不动**。
  一条一行说清"改了什么 + 为什么"，别写成设计文档；issue 号跨仓库要带包名前缀（`clusterProfiler#723`）。
- **绝不删 NEWS 条目**：精简只压缩措辞。合并小节前先 `git show HEAD:NEWS.md` 核对
  （本条已因整段误删发生过两次）。
- **结果表列变更必须追加式**：新列加尾部，原有列的顺序与集合不变，否则下游按位置取列会错位。

## 5. 已知的坑

- **`load_all()` 与 `library()` 是两个世界**：前者用工作树，后者用 `~/Rlib` 已装版本。
  验证下游（如 `clusterProfiler`）必须先 `devtools::install()`，否则测的是旧包。
- **`gseaResult` 没有 `pvalueCutoff` 槽**（只有 `enrichResult` 有；GSEA 的 cutoff 在 `@params`、
  构造时已应用）。`get_enriched()` 只服务 `enrichResult`，别写成"两类共用"。
- **`eps` 是"最小可估 p 值"的下限，不是容差**：只有进 multilevel 精算的集合碰到它，
  弱集合走 simple 分支不受影响 —— 复现 eps 问题必须用强富集集合。
- **q-value 估不出来是数据形态问题**（p 值都小且接近 → `p*n/rank` 恒定 → spline 无方差可平滑），
  只影响 `qvalue` 列，`p.adjust` 正常，所以用 `message()` 而非 `warning()`。
  ⚠ 理由**别写成"每次 print 都刷屏"**：`calculate_qvalue()` 只在**构造时**调用（7 处），
  实测 `as_enrichResult()` 发 1 条、同一对象 `print()`/`show()`/`as.data.frame()` 各 0 条。
- **NSEA 的置换只重排基因名**，`sum(geneList[common])` 是不变量，普通输入到不了退化分支；
  要测就用 `with_mocked_bindings()` 让 `.nsea_diffuse` 前几次返回 `NULL`（见 `test-nsea.R`）。
- **`dsea()` 的 design 必须是 `model.matrix()` 的数值矩阵且 rownames 与样本名一致**；
  返回的是 **data.frame**（`gsea()` 的表 + `padj` + `contrast`），不是 `gseaResult`。

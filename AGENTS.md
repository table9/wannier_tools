# AGENTS.md — 在本地 `WannierTools` 中新增 **Lindhard 函数**计算能力

> codex 多智能体协作说明文档  
> 目标：在现有本地 `WannierTools` 源码中新增 Lindhard（非相互作用极化函数/泡图）计算模块；通过输入键 `Lindhard_calc` 开启该功能，并提供网格、展宽、温度、频率等控制选项与高效并行实现。

---

## 0) 成果物清单（Definition of Done）

1. **代码层**
   - 新增模块：`src/lindhard_mod.f90`（核心实现），`src/lindhard_io.f90`（输入/输出与解析，可选）。
   - 主程序挂载：在 `driver.f90`（或现有 task 调度处）识别 `Lindhard_calc`，调用 `lindhard_mod`。
2. **输入/输出层**
   - 在 `wt.in`（或 `wannier_tools.in`）中支持以下键（示例见 §4）：
     - `Lindhard_calc = T|F`
     - `Lindhard_qmesh = NX NY NZ` 或 `Lindhard_qpath`/`Lindhard_qfile`
     - `Lindhard_mode = static|dynamic`
     - `Lindhard_omega = ω_min ω_max Nω`（dynamic 时）
     - `Lindhard_eta = 0.005  ! eV`（展宽）
     - `Lindhard_T = 0.001   ! eV（或 K，见 §4 单位约定）`
     - `Lindhard_matrix = full|CMA|diag`（矩阵元近似）
     - `Lindhard_project = none|orbital_indices`（可选投影）
     - `Lindhard_mu = auto|value` 与 `Lindhard_doping = 0.0`（可选）
   - 输出结果：
     - 静态：`chi0_static.dat`（q 与 Re/Im χ0），可选 `chi0_static.h5`
     - 动态：`chi0_w.dat` 或 `chi0_w.h5`（包含频率维度）
3. **数值正确性与性能**
   - 单元与回归测试（§6）全部通过。
   - q/k 并行（MPI/OpenMP）可用；大网格下性能线性扩展良好。
4. **文档与示例**
   - 新增 `docs/Lindhard.md`（用户手册）与 `examples/lindhard_*`（可运行示例）。
   - 本 `AGENTS.md` 作为多智能体协作的执行手册。

---

## 1) 模块目标与接口

### 1.1 计算对象
非相互作用极化函数（Lindhard）：

\[
\chi_0(\mathbf q,\omega) = \frac{1}{N_k}\sum_{\mathbf k}\sum_{n,m}
\frac{f(\varepsilon_{n\mathbf k})-f(\varepsilon_{m,\mathbf{k+q}})}
{\varepsilon_{n\mathbf k}-\varepsilon_{m,\mathbf{k+q}}+\omega+i\eta}
\,\big|\langle n,\mathbf k \,|\, \hat{M}(\mathbf q) \,|\, m,\mathbf{k+q}\rangle\big|^2
\]

- `static` 模式：\(\omega=0\)  
- `dynamic` 模式：\(\omega\) 取样  
- \(\hat{M}(\mathbf q)\) 的近似：
  - **CMA**（constant matrix element）：\(|M|^2=1\)
  - **diag**：同一 Wannier 轨道对角叠加 \(|\sum_a u^*_{n,a}(\mathbf k)\,u_{m,a}(\mathbf{k+q})|^2\)
  - **full**：包含晶格内坐标相位 \(e^{-i\mathbf q\cdot \boldsymbol\tau_a}\)：\(|\sum_a u^*_{n,a}(\mathbf k)\,u_{m,a}(\mathbf{k+q})e^{-i\mathbf q\cdot \boldsymbol\tau_a}|^2\)

### 1.2 关键输入
- TB 本征值/本征矢（从现有 `WannierTools` Hamiltonian 求解接口取得）
- \(\mu\)、T、η、k/q 网格/路径、频率网格（如启用动态）

### 1.3 关键输出
- χ0 的实部与虚部，按 q（与 ω）组织，可选投影/分辨（轨道选择时提供子块）

---

## 2) codex 多智能体编队与分工

> 统一工作分支：`feature/lindhard`  
> 统一 Issue：`#WT-CHI0`（子任务见下）

### Agent A — 规划&架构（Planner/Architect）
- **任务**：细化实现方案、接口与数据流；确认与 WT 现有模块的边界。  
- **交付**：`/docs/Lindhard.md#design` 初稿、CMake 改动清单、数据结构草图。  
- **完成标准**：得到实现同意（代码评审通过）。

### Agent B — Fortran 实现（Fortranist）
- **任务**：编写 `src/lindhard_mod.f90` 与 `src/lindhard_io.f90`；挂接主流程。  
- **重点**：数值稳定性（η 展宽、温度 Fermi 分布）、q↔k+q 索引映射、并行切分。  
- **完成标准**：示例可运行、数值测试通过。

### Agent C — 并行与性能（Perf/Parallel）
- **任务**：基于 MPI（q 层）+ OpenMP（k 层）混合并行；缓存/复用本征信息。  
- **完成标准**：规模化网格下获得近线性扩展；提供性能报告。

### Agent D — IO & 单位（I/O Wrangler）
- **任务**：输入解析、单位一致性（eV/K）、输出格式（*.dat）与兼容性。  
- **完成标准**：与旧版输入完全兼容；新增键对旧任务无副作用。

### Agent E — 测试与验证（Tester）
- **任务**：构建可验证体系与基准（§6），CI 脚本，边界条件测试。  
- **完成标准**：所有测试绿灯；误差分析记录完整。

### Agent F — 文档与示例（Docs/Examples）
- **任务**：用户手册、示例输入、可视化脚本（Python/gnuplot）。  
- **完成标准**：用户按手册 10 分钟内复现实例图。

---

## 3) 代码改动总览

### 3.1 新文件
- `src/lindhard_mod.f90`  
  - `subroutine compute_lindhard(ctx, ham, kmesh, qspec, omega_spec, chi0_out)`
  - `subroutine fill_k_eigs_cache(...)`（可选：预计算/缓存 ε、U）
  - `subroutine map_k_plus_q_index(...)`（周期边界与索引映射）
- `src/lindhard_io.f90`  
  - `subroutine read_lindhard_keywords(...)`
  - `subroutine write_chi0_ascii(...)`
- `docs/Lindhard.md`、`examples/lindhard_static_square/`、`examples/lindhard_dynamic_tb/`

### 3.2 挂载点
- `src/driver.f90`（或任务路由文件）：
  - 解析 `Lindhard_calc`
  - 当为 `T` 时调用 `compute_lindhard(...)` 并结束（或与其他任务并列）


---

## 4) 新增输入键与示例（wt.in）

> **单位约定**：  
> - 能量：**eV**（与 WT 主体一致）  
> - 温度：可用 **K** 或 **eV**，通过键区分：
>   - `Lindhard_T_unit = K|eV`（默认 K；内部换算 `k_B = 8.617333262e-5 eV/K`）

```ini
# —— 基本开关 ——
Lindhard_calc = T
Lindhard_mode = static        # static | dynamic

# —— 网格与路径（2 选 1）——
Lindhard_qmesh = 64 64 1      # 均匀 q 网格
# Lindhard_qpath = G M K G     # 或者走高对称路径（若已有路径解析器）
# Lindhard_qfile = qlist.dat   # 或外部文件指定 q 列表

# —— 频率（动态模式时需要）——
# Lindhard_omega = 0.0 0.5 251  # eV: [ω_min, ω_max, Nω]

# —— 物理参数 ——
Lindhard_eta = 0.005          # eV, 展宽
Lindhard_T = 300              # 默认单位 K
Lindhard_T_unit = K
Lindhard_mu = auto            # auto | value(eV)
Lindhard_doping = 0.0         # 每原胞电子数变化（可选）

# —— 矩阵元与投影 ——
Lindhard_matrix = full        # full | CMA | diag
Lindhard_project = none       # none | 1,2,3   (给出 Wannier 轨道索引列表)

# —— 输出控制 ——
Lindhard_out = chi0_static.dat
```

---

## 5) 核心数值流程（伪代码）

```fortran
! 在 compute_lindhard 中
call build_k_mesh(kmesh)                         ! 复用 WT 现有 k 网格
call maybe_precompute_eigs_U(kmesh, cache)       ! 可选缓存 ε_{n,k}, U_{a,n}(k)

do iq = 1, Nq   ! MPI rank 级并行
  q = qlist(iq)
  chi0 = 0.0_complex

  !$omp parallel do reduction(+:chi0)
  do ik = 1, Nk
    kpq = map_k_plus_q_index(ik, q, kmesh)      ! 周期映射到网格索引
    eig_k  = get_eigs(cache, ik)
    eig_kq = get_eigs(cache, kpq)
    U_k    = get_U(cache, ik)                   ! 本征矢 (orbitals x bands)
    U_kq   = get_U(cache, kpq)

    do n=1, nbands
      fn = fermi(eig_k(n), mu, T)
      do m=1, nbands
        fm = fermi(eig_kq(m), mu, T)
        denom = eig_k(n) - eig_kq(m) + omega + i*eta   ! static: omega=0
        if (Lindhard_matrix == 'CMA') then
          M2 = 1.0
        else
          M = sum_a( conj(U_k(a,n)) * U_kq(a,m) * phase(a,q) )  ! full/diag
          M2 = abs(M)**2
        end if
        chi0 = chi0 + (fn - fm) / denom * M2
      end do
    end do
  end do
  !$omp end parallel do

  chi0 = chi0 / Nk
  call write_out(q, chi0, outfile)
end do
```

- `phase(a,q) = exp(-i * dot(q, tau_a))`（`diag` 时置 1）
- 动态模式：对每个 `ω` 外层循环或将 `denom` 向量化
- 加速：能带裁剪（只累计靠近 μ±E_cut 的能带对），k 分块缓存

---

## 6) 验证与测试（Tester 执行）

1. **1D 最近邻 TB（半填充）**  
   - 解析特征：`q=2k_F` 处出现 Kohn 异常（cusp）  
   - 用小 η、T→0 静态 χ0(q) 重现尖峰位置  
2. **2D 方格 TB（半填充）**  
   - 验证 `full/diag/CMA` 三种矩阵元一致的极值 q 位置（幅值不同可接受）  
3. **网格收敛**  
   - 固定 q，`Nk: 40^2 → 100^2` 收敛至 1% 内  
4. **温度/展宽**  
   - 增大 T 与 η，尖峰平滑但峰位不漂移  
5. **并行一致性**  
   - 不同 MPI/OpenMP 拓扑得到数值一致（浮点误差级）  
6. **回归**  
   - 保存基准输出（*.dat/h5），后续变更自动对比

每项测试附 `examples/*` 与简易绘图脚本（Python：matplotlib；gnuplot 备选）。

---

## 7) 性能与并行建议（Perf/Parallel 执行）

- **并行切分**：MPI over q（天然解耦），OpenMP over k  
- **缓存策略**：预计算并缓存全部 k 点的 `(ε, U)`（内存允许时）；否则分块计算  
- **带裁剪**：设 `E_window`（如 0.5–1.0 eV）仅累计靠近 μ 的能带对  
- **矢量化**：对 m 循环展开；对动态模式将 `denom(ω)` 向量化  
- **I/O 降压**：默认写 ASCII 汇总

---

## 8) 风险与兼容性

- **单位混淆**：严格区分 eV/K；提供 `Lindhard_T_unit`；文档突出说明  
- **k↔k+q 映射**：必须使用规则均匀网格；若是外部 k 列表需建立哈希近邻匹配（不推荐）  
- **大内存占用**：nbands×Nk×(eig+vec) 可能很大；提供分块模式开关  
- **SOC/自旋基**：若基函数含自旋，自然并入 `U`；χ0 默认为电荷通道（总和），后续可扩展自旋通道

---

## 9) 提交与评审流程

1. `feature/lindhard`：按子任务提交原子 PR  
2. PR 模板需包含：目的、接口变更、性能数据、测试截图  
3. 必须通过：编译、示例运行、单元/回归测试、跨平台（GCC/IFX）构建

---

## 10) 开发顺序（甘特化最小闭环）

1. 解析键与主流程挂接（A/D）  
2. `CMA + static + qmesh` 最小实现（B）  
3. k↔k+q 映射/缓存与并行（C）  
4. `diag/full` 矩阵元（B）  
5. `dynamic ω` 扩展（B/C）  
6. I/O（ASCII）与示例（D/F）  
7. 基准与性能验证（E/C）  
8. 文档完善（F）

---

## 11) 用户快速上手（示例）

**静态 χ0（2D 方格 TB）**
```ini
Lindhard_calc   = T
Lindhard_mode   = static
Lindhard_qmesh  = 201 201 1
Lindhard_eta    = 0.003
Lindhard_T      = 20
Lindhard_T_unit = K
Lindhard_matrix = full
Lindhard_mu     = auto
Lindhard_out    = chi0_static.dat
```
输出：`qx qy Re(chi0) Im(chi0)`；用 `plot_chi0.py` 可直接出热力图。

**动态 χ0（一路 ω）**
```ini
Lindhard_calc   = T
Lindhard_mode   = dynamic
Lindhard_qmesh  = 81 81 1
Lindhard_eta    = 0.01
Lindhard_T      = 0.001
Lindhard_T_unit = eV
Lindhard_omega  = 0.0 0.5 251
Lindhard_matrix = diag
```

---

## 12) 附：接口与数据结构（Fortran 签名建议）

```fortran
module lindhard_types
  type LindhardInput
    logical :: do_calc
    character(len=8) :: mode      ! 'static'/'dynamic'
    integer :: qmesh(3)
    real(8) :: eta, temp, mu
    character(len=8) :: T_unit    ! 'K'/'eV'
    character(len=8) :: matrix    ! 'CMA'/'diag'/'full'
    real(8) :: omega_min, omega_max
    integer :: nomega
    integer, allocatable :: proj_list(:)
  end type
end module

module lindhard_mod
contains
  subroutine compute_lindhard(ctx, ham, lind_in)
  ! ctx/ham 复用 WT 的上下文与哈密顿接口
  end subroutine
end module
```

---

## 13) 维护与扩展路线

- 追加 **自旋/电荷通道分解** 与 **轨道分辨 χ0**  
- 支持 **四面体积分**（减少 η 依赖）  
- 对接 **RPA**：\(\chi = \chi_0(1-U\chi_0)^{-1}\)（简单 Hubbard-U 近似）

---

### 提醒
- 本改动对原有 `WannierTools` 任务保持默认兼容：`Lindhard_calc = F` 时不触发任何新逻辑。  

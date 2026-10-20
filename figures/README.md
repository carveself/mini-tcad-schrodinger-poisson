# 默认计算数据

8 张 PNG、10 个 CSV、1 个 NPZ、3 个 JSON，均来自理想模型计算，无实测数据。

`bias_sweep.csv` 的经典支路没有离散子带，相关列为空；量子波函数在无限壁为零。`quantum_state.npz` 保存 0.8 V 的 z（m）、电势（V）、密度（m⁻³）、能级（eV）与波函数（m⁻¹/²）。

面电荷单位 cm⁻²；电容为两栅总电荷对共同栅压的导数，单位 F/m²。最近界面质心是 `积分 n*min(z,tSi-z)/积分 n`，不是全膜平均 z。量子电容的冻结谱定义见主 README，不能直接当成栅端电容。

`state_truncation.csv` 在冻结电势下比较保留态数；`mixing_check.csv` 和 `iteration_traces.json` 检查真实电荷残差。谱的 FEM 参考只检查冻结势阱的能级。默认静态报告、参数与适用范围见 `docs/`；改变参数后保留自己的记录。

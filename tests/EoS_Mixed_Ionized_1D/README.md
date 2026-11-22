# Mixed / Ionized EoS 回归算例

三个 1D 算例共用一个初始化文件，`HD_NCOMP=2`，使用 CGS 单位。这里不调用 Saha，电离态作为独立组分输运。

| 参数文件 | 初始状态 | 检查内容 |
| --- | --- | --- |
| `parameters_mixed_contact.txt` | H I / He I，周期平滑密度与组分波，P=100、v=1e5，演化 200 s | Mixed 的 mu、温度、零自由电子数；接触输运、压力平衡、质量和能量守恒 |
| `parameters_ionized_contact.txt` | H I / H II，同样的周期接触输运 | 自由电子对 mu 的贡献、电离能随 H II 输运、压力平衡及组分守恒 |
| `parameters_ionized_heating.txt` | rho=1e-10，T=1e4，X(H I)=0.65，X(H II)=0.35，均匀加热 0.2 s | 两套能量同步增加 rho*q*dt；压力和温度变化；源项时间尺度使用含电离能的内能 |

直接运行一个算例：

```bash
cd tests/EoS_Mixed_Ionized_1D
make -j4
./flare -par parameters_mixed_contact.txt
./flare -par parameters_ionized_contact.txt
./flare -par parameters_ionized_heating.txt
```

自动运行并校验：

```bash
python3 check.py --exe ./flare
```

脚本为两个接触算例分别运行 RK2、LTS RK2；为加热算例运行 Euler、RK2、LTS Euler、LTS RK2、Implicit，共 9 组。运行时生成的参数和输出默认保存在临时目录，末尾打印路径；也可以用 `--output-root /tmp/my-eos-check` 指定目录。

检查包括：终止时间、质量/各组分/能量积分、组分和、恒定接触压力与速度、输运后的组分波形，以及加热的解析能量和压力。接触算例使用 PCM，组分波存在预期的数值扩散，L1 误差要求小于 0.025。

用户输出钩子独立计算粒子数和能量关系，检查当前缓存中的温度、mu、电子数密度，并把最后一次快照的这些量写入 `thermodynamics_last.csv`，便于核对。框架原有 CSV 格式没有改变。

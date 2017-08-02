# CPC2017_WENO
WENO格式（Weighted Essentially Non-Oscillatory schemes) 是基于ENO格式（Essentially Non-Oscillatory schemes）发展而来的一种求解双曲守恒律组的高精度高分辨率广义Godunov格式，适合于求解高密度比流体界面不稳定等具有强间断又具有大面复杂流动结构的问题。本程序由中科院力学所李新亮老师团队开发的OpenCFD程序精简而来，在此对李老师团队致以诚挚的谢意。 程序包含4个源文件：main.f90, Interfaces.f90， Weno.f90，parameters.f90。 其中main.f90是程序主要框架；Weno.f90为程序主要算法实现；Interfaces.f90 中为通信函数和IO函数的封装； parameters.f90中包含了程序部分常量声明。

# CPC2017_WENO
WENO格式（Weighted Essentially Non-Oscillatory schemes) 是基于ENO格式（Essentially Non-Oscillatory schemes）发展而来的一种求解双曲守恒律组的高精度高分辨率广义Godunov格式，适合于求解高密度比流体界面不稳定等具有强间断又具有大面复杂流动结构的问题。本程序由中科院力学所李新亮老师团队开发的OpenCFD程序精简而来，在此对李老师团队致以诚挚的谢意。 程序包含4个源文件：main.f90, Interfaces.f90， Weno.f90，parameters.f90。 其中main.f90是程序主要框架；Weno.f90为程序主要算法实现；Interfaces.f90 中为通信函数和IO函数的封装； parameters.f90中包含了程序部分常量声明。
# 编译：
直接在code目录下 make即可。  

$make  

mpif90  -cpp -DTEST -DSWAP -DCHECK -c-OPT:ieee_arith=1 -c parameters.f90  

mpif90  -cpp -DTEST -DSWAP -DCHECK -c-OPT:ieee_arith=1 -c main.f90  

mpif90  -cpp -DTEST -DSWAP -DCHECK -c -OPT:ieee_arith=1-c Weno.f90  

mpif90  -cpp -DTEST -DSWAP -DCHECK -c-OPT:ieee_arith=1 -c Interfaces.f90  

mpif90  -O3  -o../run/weno7.out parameters.o main.o Weno.o Interfaces.o  

# 运行：
本次比赛设置有三个不同的算例，分别为exp1、exp2和exp3。每个目录下均有./run.sh脚本，直接运行该脚本可以将编译成功的可执行文件拷贝到当前目录下，并提交程序到队列中。
脚本内容如下：
bsub -o run.log -b -q  q_sw_expr -n 4 -np 4 -cgsp 64 -share_size 6000 -host_stack 512 -priv_size 16  -pe_stack 3 ./weno7.out

[初赛题目介绍](http://mp.weixin.qq.com/s?__biz=MzU5MzAzMDM4Nw==&mid=2247483709&idx=1&sn=7dd0492250bc02a821fed738b57bc6e5&chksm=fe17fca1c96075b7e922923e0248db0d57dfcad4225104f8b5e3ff5afd9c6fc66f6097e50a99&mpshare=1&scene=23&srcid=0801ObC74KgWYwVL9cSF1dPe#rd)  
[代码与算例更新](http://mp.weixin.qq.com/s?__biz=MzU5MzAzMDM4Nw==&mid=2247483713&idx=1&sn=4a6b31d0ce994b8f46d2e207fca72264&chksm=fe17fcddc96075cb3970d75aa10531cb504b6f9e168dd7dda75531fbf45ef593ca9909ff0a2b&mpshare=1&scene=23&srcid=0801Y2rKyLftFrIExibMvE9t#rd)  
[校验程序](http://mp.weixin.qq.com/s?__biz=MzU5MzAzMDM4Nw==&mid=2247483713&idx=2&sn=a6bace7ea0d3af5e4ad151b8f15beb53&chksm=fe17fcddc96075cbddc7f7f3b2ad89a5e2cb26484a38166438ccd7e81a254c0072bc87080b60&mpshare=1&scene=23&srcid=0801FRRBDV29vKNDQ5lOt99G#rd)  

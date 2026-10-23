#pragma once
#include <vector>
#include <string>
#include "hd_quantities.h"
//Global Domain
//1-D resize(N)
// 所有 .cpp 共用同一个数组（需要 C++17），初始为空。
// GridCells.resize(N)：调整单元数，新增单元调用 GridValue 默认构造函数。
// GridCells[i]：访问第 i 个单元，下标从 0 开始。
inline std::vector<Grid> Domain;
inline int ndim,ncon,ncomp,ncon_custorm,nvar;  //ndim - dimension     ncon - built-in conservative params number, it varies with dimension and mhd involvement. ncon_custom is user-defined conservatives ncomp is the composition
inline int nx,ny;   //active domain + ghostcells (nxa + 2*nxg,...)
inline int nxg,nyg; //ghost cell width on EACH side
inline int nxa,nya; //active domain
inline double tmax,dt; //currently fixed, should be adaptive
inline int tout;    //how many timesteps for an output file?
inline double Xa,dxa; //Xa is the total length of the active region; dxa = Xa/nxa
inline std::string boundary_conditions; //implement it later
inline int output_every;    //after how many steps we write an output
inline std::string output_name; // name of output prefix

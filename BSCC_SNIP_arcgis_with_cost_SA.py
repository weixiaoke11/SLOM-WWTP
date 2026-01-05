# coding=GBK
# noinspection PyBroadException
# Copyright 2014  Swiss Federal Institute of Aquatic Science and Technology
#
# This file is part of SNIP (Sustainable Network Infrastructure Planning)
# SNIP is used for determining the optimal degree of centralization for waste
# water infrastructures. You find detailed information about SNIP in Eggimann et al. (2014).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#    
# The Djikstra and a* algorithm are adapted from Hetland (2010).
# The algorithm is developed for Python 2.7 and ArcGIS 10.2
#
# Literature
# ==========
# Eggimann Sven, Truffer Bernhard, Maurer Max (2015): To connect or not to connect? 
# Modelling the optimal degree of centralisation for wastewater infrastructures.   
# Water Research, 84, 218-231. Link: http://www.sciencedirect.com/science/article/pii/S004313541530107X
#
# Hetland M.L. (2010): Python Algorithms. Mastering Basic Algorithms in the Python Language. apress.
#
# Contact:   sven.eggimann@eawag.ch
# Version    1.0
# Date:      1.07.2015
# Autor:     Eggimann Sven


import gc            # Don't allow garbage collection
import pickle
import os
import numpy
import sys
import argparse
from BSCC_SNIP_functions_with_cost_SA import *                                # Import Functions


if __name__ == '__main__':
    gc.disable()
    parser = argparse.ArgumentParser()

    # 实验序号，含shp名
    parser.add_argument("--exp", type=str, default="1")
    # 城市文件，含shp名
    parser.add_argument("--city_shp", type=str, default="1.shp")
    # 城市保存的文件夹序号
    parser.add_argument("--city_folder", type=str, default="1")
    # 城市文名称，检索参数用
    parser.add_argument("--city_unique", type=str, default="北京市")

    # 扩展模块初始温度与降温速度
    parser.add_argument("--sa_t_ori_e", type=float, default=10000)
    parser.add_argument("--sa_alpha_ori_e", type=float, default=0.95)
    # 合并模块初始温度与降温速度
    parser.add_argument("--sa_t_ori_m", type=float, default=10000)
    parser.add_argument("--sa_alpha_ori_m", type=float, default=0.95)
    # 合并模块的温度是否与扩展共同减小
    parser.add_argument("--sa_co_vary", type=bool, default=False)

    # 需达到Z的大小
    parser.add_argument("--z_to_reach", type=float, default=1)

    args = parser.parse_args()

    exp = args.exp
    city_shp = args.city_shp
    city_folder = args.city_folder
    city_unique = args.city_unique

    sa_t_ori_e = args.sa_t_ori_e
    sa_alpha_ori_e = args.sa_alpha_ori_e
    sa_t_ori_m = args.sa_t_ori_m
    sa_alpha_ori_m = args.sa_alpha_ori_m
    sa_co_vary = args.sa_co_vary

    z_to_reach = args.z_to_reach

    print("city:{}-{}, folder:{}".format(city_unique, city_shp, city_folder))
    out_data = f"./out_data_SA_{exp}"
    # 生成文件夹
    if not os.path.exists(out_data):
        os.mkdir(out_data)
    # 生成文件夹
    if not os.path.exists("{}/{}".format(out_data, city_shp[0:-4])):
        os.mkdir("{}/{}".format(out_data, city_shp[0:-4]))
    # 生成进度条文件夹
    if not os.path.exists(r"./log"):
        os.mkdir(r"./log")

    outListFolder = "{}/{}".format(out_data, city_shp[0:-4]) + "/"

    (_, forSNIP, streetGraph, startnode, edgeList, streetVertices, rasterSize, buildPoints, buildings, rasterPoints,
     InputParameter, aggregatetPoints, StartXY) = pickle.load(open(r"./in_data/{}/data.dat".format(city_shp[0:-4]),"rb"))

    InputParameter["parallel"] = False

    # Run SNIP
    (ExpansionTime, MergeTime, sewers, pointsPrim, WWTPs, wtpstodraw, pumpList, edgeList,
     completePumpCosts, completeWWTPCosts, completePublicPipeCosts, totalSystemCosts,
     buildings, buildPoints, aggregatetPoints, _) = SNIP(
        0, outListFolder, 1, forSNIP, 1, streetGraph, startnode,edgeList, streetVertices, rasterSize, buildPoints,
        buildings, rasterPoints, InputParameter, aggregatetPoints, z_to_reach=z_to_reach,
        sa_t_ori_e=sa_t_ori_e, sa_alpha_ori_e=sa_alpha_ori_e,
        sa_t_ori_m=sa_t_ori_m, sa_alpha_ori_m=sa_alpha_ori_m, sa_co_vary=sa_co_vary,)

    pickle.dump((ExpansionTime, MergeTime, sewers, pointsPrim, WWTPs, wtpstodraw, pumpList, edgeList, completePumpCosts,
                 completeWWTPCosts, completePublicPipeCosts, totalSystemCosts, buildings, buildPoints, aggregatetPoints),
                open("{}/{}/data_out.dat".format(out_data, city_shp[0:-4]), "wb"))


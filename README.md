# part1

### process_group_and_table函数：

**功能**：匹配两个表（丰度表和分组信息表），两个表共有的保留，单独的不要

### filter_OTU2函数

**功能**：根据流行度过滤

### filter_OTU函数

**功能**：根据top值来筛选OTU

### corMicro函数

**功能**：构建相关性矩阵

**参数**：1.丰度表table1	2.相关性阈值r.threshold/p.thresholde	3.构建相关性的方法method

### node_properties函数

**功能**：计算网络的全局和局部属性

**参数**：1.图对象g，格式为graph	2.网络聚类方法clu_method

### plotnetwork可视化网络

**参数**：1.图对象g，格式为graph	2.网络聚类方法clu_method	3.tag=sel_group表示分组标签	4.nod_cluster网络低于某值为灰色	5.output为输出路径

### global_pro_compare函数

**功能**：和随机网络比较

**参数**：1.图对象g，格式为graph	2.step随机生成次数	3.netName网络分组标签	4.ncpus所用cpu核数

### Feature.Random.removal函数

**功能**：模拟移除相应物种后的鲁棒性、复杂性、脆弱性、稳定性

# part2

### process_group_and_table函数：

**功能**：匹配两个表（丰度表和分组信息表），两个表共有的保留，单独的不要

### filter_OTU2函数

**功能**：根据流行度过滤

### filter_OTU函数

**功能**：根据top值来筛选OTU

### merge_bio函数

**功能**：合并两个表格数据

**参数**：1.表1table1	2.表2table2	3.table1表命名	4.table2表命名

### corMicro函数

**功能**：构建相关性矩阵

**参数**：1.丰度表table1	2.相关性阈值r.threshold/p.thresholde	3.构建相关性的方法method

### node_properties函数

**功能**：计算网络的全局和局部属性

**参数**：1.图对象g，格式为graph	2.网络聚类方法clu_method

### bridge_network桥接网络可视化

**参数**：1.混合的丰度矩阵	2.混合的相关性矩阵	3.OTU1分组命名

# part3

### process_group_and_table函数：

**功能**：匹配两个表（丰度表和分组信息表），两个表共有的保留，单独的不要

### filter_OTU2函数

**功能**：根据流行度过滤

### sspcc_cal函数

**功能**：将所有（n-1）样品构建为背景网络，进行单样品网络构建

**参数**：1.丰度表

### sspcc_cal3函数

**功能**：根据分组信息将每组内的（n-1）样品设为背景网络

**参数**：1.丰度表	2.分组信息	3.分组列名vscol1

###  sspcc_cal2函数

**功能**：选取某一组为对照组构建单样品网络

！！！**参数**：1.丰度表	2.分组信息	3.分组列名vscol1	4.对照组名

### make_one_edgelist_file_general函数：

**功能**：合并矩阵

### scale_network函数：

**功能**：标准化网络矩阵

**参数**：1.相关性矩阵 	2.是否保存

### calculate_sum_of_weight函数

**功能**：计算单样品网络的节点权重

**参数**：过滤后的单样品网络矩阵

### PCA_draw函数：

**功能**：计算PCA/PCOA的可视化

**参数**：1.SOW单样品网络点的权重矩阵	2.cor_top单样品网络边的权重矩阵	3.vscol2分组列名2	4.ShowType显示PCA可视化参数，如centroid	5.mode选择PCA和PCOA	6.offset是否合并第一列和第二列信息

### Differential_nodes函数：

**功能**：计算单样品网络差异节点和差异边以及可视化

**参数**：是否可视化

### determineCharacteristics函数：

**功能**：计算单样品网络属性

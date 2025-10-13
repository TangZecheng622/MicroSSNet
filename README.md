# MicroSSNet: An R package for microbial network construction and analysis at the single-sample and aggregated levels.
Network analysis is a fundamental tool for elucidating microbial interactions, which are crucial for understanding the mechanisms that shape ecosystem structure and function. Traditional aggregated (co-abundance/co-occurrence) network approaches that infer pairwise relationships among biological entities from large sample collections often overlook sample-specific interaction patterns. To address this, we developed MicroSSNet, an R package designed for analyzing microbial networks at both the aggregated and single-sample levels. These results highlight the importance of sample-specific interaction patterns and demonstrate that SSN-based approaches provide complementary insights into aggregated network and abundance-based approaches. MicroSSNet offers a robust framework for constructing and analyzing single-sample microbial networks, advancing microbiome research at both the individual and community scales. The package is freely available on GitHub (https://github.com/TangZecheng622/MicroSSNet).

# Examples of visualizations.
<img width="729" height="1050" alt="image" src="https://github.com/user-attachments/assets/eb9f4c5f-c334-4999-9b6c-3ffd07a5516e" />

## 1.aggregated network
<img width="1095" height="779" alt="image" src="https://github.com/user-attachments/assets/f6caad84-6ff7-4b61-b23d-cce7570a52a8" />

## 2.single-sample network
<img width="1107" height="1036" alt="image" src="https://github.com/user-attachments/assets/166b0787-0140-4fa7-9a9f-8286de89a747" />


# Installation

~~~
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("limma")
install("remotes")
remotes::install_github("zdk123/SpiecEasi")
remotes::install_github("TangZecheng622/MicroSSNet")
~~~



# aggregation_netpipeline

~~~R
aggregated_network = aggregation_netpipeline (
    table1,
    tax = NULL,
    cor_table_list = NULL,
    group_df = NULL,
    vscol1 = "Group",
    pvl.threshold = 0.5,
    top = NULL,
    r.threshold = 0.8,
    p.threshold = 0.05,
    method = "pearson",
    R = 10,
    ncpus = 4,
    clu_method = "cluster_fast_greedy",
    p.adj = "fdr",
    step = 10,
    node.cluster = 5,
    rm.p.list = seq(0, 0.5, by = 0.1),
    zipi = TRUE,
    compare_net = TRUE,
    calculate_vul = TRUE,
    calculate_rob = TRUE,
    calculate_cpx = TRUE,
    output_dir = getwd())
~~~

### Mainly use functions:



### process_group_and_table function：

**Function**: Align two tables (an abundance table and a group information table) by matching their shared entries. Only the common entries present in both tables are retained; any entries unique to either table are discarded.

### filter_OTU2 function

**Function**: Perform filtering by retaining entries according to their prevalence.

### filter_OTU function

**Function**: Select OTUs by sorting them based on their top values and retaining the top entries.

### corMicro function

**Function**: Construct a correlation matrix.  
 **Parameters**:

1. Abundance table (`table1`)
2. Correlation thresholds (`r.threshold` / `p.threshold`)
3. Method used for constructing the correlation matrix (`method`)

### node_properties function

**Function**: Calculate global and local properties of the network.  
 **Parameters**:

1. Graph object (`g`), in `graph` format
2. Network clustering method (`clu_method`)

### plotnetwork Visualize the network

**Parameters**:

1. Graph object (`g`), in `graph` format
2. Network clustering method (`clu_method`)
3. `tag` = `sel_group`, indicating the grouping label
4. `nod_cluster`, threshold below which nodes are colored gray
5. `output`, path for saving the output

### global_pro_compare function

**Function**: Compare with random networks.  
 **Parameters**:

1. Graph object (`g`), in `graph` format
2. `step`, A randomly generated number indicating the number of times to generate random networks.
3. `netName`, network group label
4. `ncpus`, number of CPU cores used

### Feature.Random.removal function

**Function**: Simulate and evaluate the resulting robustness, complexity, vulnerability, and stability of the network following the removal of specific species.

# bipartite_netpipeline

~~~R
single_sample_network = bipartite_netpipeline (
    table1,
    table2,
    table1name = "domain1",
    table2name = "domain2",
    tax = NULL,
    cor_table_list = NULL,
    group_df = NULL,
    vscol1 = "group",
    pvl.threshold = 0.5,
    top = NULL,
    r.threshold = 0.8,
    p.threshold = 0.05,
    method = "pearson",
    R = 10,
    ncpus = 4,
    clu_method = "cluster_fast_greedy",
    p.adj = "BH",
    step = 10,
    rm.p.list = seq(0, 0.5, by = 0.1),
    zipi = TRUE,
    compare_net = TRUE,
    calculate_vul = TRUE,
    calculate_rob = TRUE,
    calculate_cpx = TRUE,
    output_dir = getwd())
~~~

### Mainly use functions:


### process_group_and_table function：

**Function**: Retain only shared entries between an abundance table and a grouping information table.

### filter_OTU2 function

**Function**: Perform filtering by retaining entries according to their prevalence.

### filter_OTU function

**Function**: Select OTUs by sorting them based on their top values and retaining the top entries.

### merge_bio function

**Function**: Merge two data tables.  
 **Parameters**:

1. `table1`, the first table
2. `table2`, the second table
3. `table1` name assignment
4. `table2` name assignment

### corMicro function

**Function**: Construct a correlation matrix.  
 **Parameters**:

1. Abundance table (`table1`)
2. Correlation thresholds (`r.threshold` / `p.threshold`)
3. Method for constructing the correlation matrix (`method`)

### node_properties function

**Function**: Calculate global and local properties of the network.  
 **Parameters**:

1. Graph object (`g`), in `graph` format
2. Network clustering method (`clu_method`)

### bridge_network Visualization

**Parameters**:

1. Mixed abundance matrix (`abundance matrix`)
2. Mixed correlation matrix (`correlation matrix`)
3. Assigned group name for `OTU1`

# ssn_pipeline

~~~R
tab_s_pip =ssn_pipeline (
    table1,
    group_df = NULL,
    vscol1,
    vscol2 = vscol1,
    ssn_method = "ssPCC",
    control = "pergroup",
    log = TRUE,
    top = NULL,
    pvl_threshold = 0.5,
    r_threshold = 0,
    scale = TRUE,
    save = TRUE,
    pca = TRUE,
    pcoa = TRUE,
    dis_method = "euclidean",
    binary = FALSE,
    showType = NULL,
    limma = TRUE,
    property = TRUE
) 
~~~

### Mainly use functions:


### process_group_and_table function：

**Function**: Match two tables (an abundance table and a grouping information table), retaining only the shared entries and discarding those unique to either table.

### filter_OTU2 function

**Function**: Perform filtering based on prevalence.

### sspcc_cal function

**Function**: Build a background network based on all (n-1) samples, and generate a single sample network.
 **Parameters**:

Abundance table

### sspcc_cal3 function

**Function**: Based on the grouping information, construct a background network for each group using (n-1) samples.
 **Parameters**:

1. Abundance table
2. Grouping information table
3. Grouping column name (`vscol1`)

###  sspcc_cal2 function

**Function**: Construct single-sample networks by selecting a specified group as the control group.
 **Parameters**:

1. Abundance table
2. Grouping information table
3. Name of the grouping column (`vscol1`)
4. Name of the control group

### make_one_edgelist_file_general function：

**Function**: Merge multiple matrices into a single matrix.

### scale_network function：

**Function**: Perform normalization of the network correlation matrix.
 **Parameters**:

1. Correlation matrix
2. Option to save the normalized matrix

### calculate_sum_of_weight function

**Function**: Calculate the node weights of single-sample networks based on the filtered network matrix.
 **Parameters**:

Filtered single-sample network matrix

### PCA_draw function：

**Function**: Perform and visualize PCA/PCoA analysis.
 **Parameters**:

1. `SOW`, node weight matrix of the single-sample network
2. `cor_top`, edge weight matrix of the single-sample network
3. `vscol2`, name of the second grouping column
4. `ShowType`, parameters for PCA visualization (e.g., centroid)
5. `mode`, selection between PCA and PCoA
6. `offset`, whether to merge the first and second column information

### Differential_nodes function：

**Function**: Calculate differential nodes and edges in single-sample networks and visualize the results.
 **Parameters**:

Visualization option (`True` or `False`)

### determineCharacteristics function：

**Function**: Calculate the network properties of single-sample networks.


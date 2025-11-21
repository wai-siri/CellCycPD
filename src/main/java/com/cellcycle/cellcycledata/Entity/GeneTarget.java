package com.cellcycle.cellcycledata.Entity;

import com.baomidou.mybatisplus.annotation.TableName;
import lombok.Data;

@Data
@TableName("gene_target")
public class GeneTarget {
    private String gene;
    private String stage;
    private String drug;
    private String cancer;
    private String pmid;
}

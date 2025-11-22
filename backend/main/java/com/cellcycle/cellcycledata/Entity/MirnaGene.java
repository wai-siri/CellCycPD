package com.cellcycle.cellcycledata.Entity;

import com.baomidou.mybatisplus.annotation.TableName;
import lombok.Data;

/**
 * MirnaGene实体类，对应数据库中的mirna_gene表
 */
@Data
@TableName("mirna_gene")
public class MirnaGene {
    // 所有字段都使用与数据库列名完全一致的命名
    private String miRNA;
    private String gene;
}

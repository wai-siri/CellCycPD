package com.cellcycle.cellcycledata.Entity;

import com.baomidou.mybatisplus.annotation.TableName;
import lombok.Data;

@Data
@TableName("ppi")
public class PPIData {
    // 所有字段都使用与数据库列名完全一致的命名
    private String node1;
    private String node2;
    private String node1_string_id;
    private String node2_string_id;
    private Integer neighborhood_on_chromosome;
    private Double gene_fusion;
    private Double phylogenetic_cooccurrence;
    private Double homology;
    private Double coexpression;
    private Double experimentally_determined_interaction;
    private Double database_annotated;
    private Double automated_textmining;
    private Double combined_score;
}

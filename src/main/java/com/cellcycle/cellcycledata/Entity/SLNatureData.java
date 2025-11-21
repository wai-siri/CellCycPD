package com.cellcycle.cellcycledata.Entity;

import com.baomidou.mybatisplus.annotation.TableField;
import com.baomidou.mybatisplus.annotation.TableName;
import lombok.Data;

/**
 * 合成致死自然数据实体类
 * 对应数据库中的sl_nature表
 */
@Data
@TableName("sl_nature")
public class SLNatureData {
    
    /**
     * 基因A
     */
    @TableField("Gene_A")
    private String geneA;
    
    /**
     * 基因A的阶段
     */
    @TableField("stage_Gene_A")
    private String stageGeneA;
    
    /**
     * 基因B
     */
    @TableField("Gene_B")
    private String geneB;
    
    /**
     * 基因B的阶段
     */
    @TableField("stage_Gene_B")
    private String stageGeneB;
    
    /**
     * GEMINI敏感性
     */
    @TableField("GEMINI_sensitive")
    private String geminiSensitive;
    
    /**
     * 细胞系
     */
    @TableField("Cell_line")
    private String cellLine;
}

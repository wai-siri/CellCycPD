package com.cellcycle.cellcycledata.Entity;

import com.baomidou.mybatisplus.annotation.TableName;
import lombok.Data;

/**
 * 合成致死(Synthetic Lethality)数据实体类
 * 对应数据库中的sl表
 */
@Data
@TableName("sl")
public class SLData {
    
    /**
     * 基因A
     */
    private String geneA;
    
    /**
     * 基因B
     */
    private String geneB;
}

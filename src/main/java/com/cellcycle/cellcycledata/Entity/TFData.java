package com.cellcycle.cellcycledata.Entity;

import com.baomidou.mybatisplus.annotation.TableName;
import lombok.Data;

/**
 * 转录因子-基因互作数据实体类
 * 对应数据库中的tf表
 */
@Data
@TableName("tf")
public class TFData {
    // 所有字段都使用与数据库列名完全一致的命名
    private String tf;
    private String target_gene;
    private String interaction;
}

package com.cellcycle.cellcycledata.Entity;

import com.baomidou.mybatisplus.annotation.TableName;
import lombok.Data;

@Data
@TableName("drug_information")
public class DrugInformation {
    private String drug;
    private String information;
}

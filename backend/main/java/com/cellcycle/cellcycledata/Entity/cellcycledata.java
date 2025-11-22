package com.cellcycle.cellcycledata.Entity;

import com.baomidou.mybatisplus.annotation.TableName;
import lombok.Data;

@Data
@TableName("cellcyclegenedata")
public class cellcycledata {
    private String Gene_Symbol;
    private String Description;
    private String Category;
    private String Species;
    private String Gene_ID;
    private String Related_Gene_Symbol;
    private String UniProt_ID;
    private String Article_Title;
    private String PMID;
    private String Year_of_Publication;
    private String Method;
    private String stage;
  
    public String getGene_Symbol(){return Gene_Symbol;}
}
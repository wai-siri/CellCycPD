package com.cellcycle.cellcycledata.Dao;

import org.apache.ibatis.annotations.Mapper;
import com.cellcycle.cellcycledata.Entity.DrugInformation;

import java.util.List;

@Mapper
public interface DrugInformationDao {
    // 查询所有药物信息
    List<DrugInformation> selectAllDrugInformation();
    
    // 根据药物名称查询信息
    DrugInformation selectByDrug(String drug);
}

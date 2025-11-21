package com.cellcycle.cellcycledata.Dao;

import org.apache.ibatis.annotations.Mapper;
import com.cellcycle.cellcycledata.Entity.GeneTarget;

import java.util.List;

@Mapper
public interface GeneTargetDao {
    // 查询所有基因靶点数据
    List<GeneTarget> selectAllGeneTargets();
    
    // 根据基因名查询
    List<GeneTarget> selectByGene(String gene);
    
    // 根据阶段查询
    List<GeneTarget> selectByStage(String stage);
}

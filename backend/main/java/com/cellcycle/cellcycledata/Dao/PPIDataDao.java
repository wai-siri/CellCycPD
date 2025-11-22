package com.cellcycle.cellcycledata.Dao;

import org.apache.ibatis.annotations.Mapper;
import org.apache.ibatis.annotations.Param;
import org.apache.ibatis.annotations.Select;
import com.cellcycle.cellcycledata.Entity.PPIData;

import java.util.List;
import java.util.Map;

@Mapper
public interface PPIDataDao {
    /**
     * 查询所有蛋白质互作数据
     */
    @Select("SELECT * FROM ppi")
    List<PPIData> selectAllPPIData();
    
    /**
     * 调试：查询表结构
     */
    @Select("SHOW COLUMNS FROM ppi")
    List<Map<String, Object>> showColumns();
    
    /**
     * 调试：原始SQL查询
     */
    @Select("SELECT node1, node2, node1_string_id, node2_string_id, combined_score FROM ppi LIMIT 1")
    List<Map<String, Object>> selectRawData();
    
    /**
     * 根据最低分数阈值查询蛋白质互作数据
     */
    @Select("SELECT * FROM ppi WHERE combined_score >= #{minScore}")
    List<PPIData> selectPPIDataByMinScore(@Param("minScore") Double minScore);
    
    /**
     * 查询与特定基因相关的蛋白质互作数据
     */
    @Select("SELECT * FROM ppi WHERE node1 = #{gene} OR node2 = #{gene}")
    List<PPIData> selectPPIDataByGene(@Param("gene") String gene);
    
    /**
     * 根据基因和最低分数阈值查询蛋白质互作数据
     */
    @Select("SELECT * FROM ppi WHERE (node1 = #{gene} OR node2 = #{gene}) AND combined_score >= #{minScore}")
    List<PPIData> selectPPIDataByGeneAndMinScore(@Param("gene") String gene, @Param("minScore") Double minScore);
}

package com.cellcycle.cellcycledata.Dao;

import com.cellcycle.cellcycledata.Entity.SLNatureData;
import org.apache.ibatis.annotations.Mapper;
import org.apache.ibatis.annotations.Param;
import org.apache.ibatis.annotations.Select;

import java.util.List;
import java.util.Map;

/**
 * 合成致死自然数据(SL Nature)数据访问接口
 */
@Mapper
public interface SLNatureDataDao {
    
    /**
     * 获取所有合成致死自然数据
     * @return 所有合成致死自然数据列表
     */
    @Select("SELECT Gene_A AS geneA, stage_Gene_A AS stageGeneA, Gene_B AS geneB, stage_Gene_B AS stageGeneB, GEMINI_sensitive AS geminiSensitive, Cell_line AS cellLine FROM sl_nature")
    List<SLNatureData> selectAllSLNatureData();
    
    /**
     * 调试：查询表结构
     */
    @Select("SHOW COLUMNS FROM sl_nature")
    List<Map<String, Object>> showColumns();
    
    /**
     * 调试：原始SQL查询
     */
    @Select("SELECT Gene_A AS geneA, stage_Gene_A AS stageGeneA, Gene_B AS geneB, stage_Gene_B AS stageGeneB, GEMINI_sensitive AS geminiSensitive, Cell_line AS cellLine FROM sl_nature LIMIT 1")
    List<Map<String, Object>> selectRawData();
    
    /**
     * 根据基因名称查询相关的合成致死自然数据
     * @param gene 基因名称
     * @return 与指定基因相关的合成致死自然数据列表
     */
    @Select("SELECT Gene_A AS geneA, stage_Gene_A AS stageGeneA, Gene_B AS geneB, stage_Gene_B AS stageGeneB, GEMINI_sensitive AS geminiSensitive, Cell_line AS cellLine FROM sl_nature WHERE Gene_A = #{gene} OR Gene_B = #{gene}")
    List<SLNatureData> selectSLNatureDataByGene(@Param("gene") String gene);
}

package com.cellcycle.cellcycledata.Dao;

import com.cellcycle.cellcycledata.Entity.SLData;
import org.apache.ibatis.annotations.Mapper;
import org.apache.ibatis.annotations.Param;
import org.apache.ibatis.annotations.Select;

import java.util.List;
import java.util.Map;

/**
 * 合成致死(Synthetic Lethality)数据访问接口
 */
@Mapper
public interface SLDataDao {
    
    /**
     * 获取所有合成致死数据
     * @return 所有合成致死数据列表
     */
    @Select("SELECT geneA, geneB FROM sl")
    List<SLData> selectAllSLData();
    
    /**
     * 调试：查询表结构
     */
    @Select("SHOW COLUMNS FROM sl")
    List<Map<String, Object>> showColumns();
    
    /**
     * 调试：原始SQL查询
     */
    @Select("SELECT geneA, geneB FROM sl LIMIT 1")
    List<Map<String, Object>> selectRawData();
    
    /**
     * 根据基因名称查询相关的合成致死数据
     * @param gene 基因名称
     * @return 与指定基因相关的合成致死数据列表
     */
    @Select("SELECT geneA, geneB FROM sl WHERE geneA = #{gene} OR geneB = #{gene}")
    List<SLData> selectSLDataByGene(@Param("gene") String gene);
}

package com.cellcycle.cellcycledata.Dao;

import org.apache.ibatis.annotations.Mapper;
import org.apache.ibatis.annotations.Param;
import org.apache.ibatis.annotations.Select;
import com.cellcycle.cellcycledata.Entity.TFData;

import java.util.List;
import java.util.Map;

/**
 * 转录因子-基因互作数据访问接口
 */
@Mapper
public interface TFDataDao {
    /**
     * 查询所有转录因子-基因互作数据
     */
    @Select("SELECT * FROM tf")
    List<TFData> selectAllTFData();
    
    /**
     * 调试：查询表结构
     */
    @Select("SHOW COLUMNS FROM tf")
    List<Map<String, Object>> showColumns();
    
    /**
     * 调试：原始SQL查询
     */
    @Select("SELECT tf, target_gene, interaction FROM tf LIMIT 1")
    List<Map<String, Object>> selectRawData();
    
    /**
     * 查询与特定转录因子相关的互作数据
     */
    @Select("SELECT * FROM tf WHERE tf = #{tf}")
    List<TFData> selectTFDataByTF(@Param("tf") String tf);
    
    /**
     * 查询与特定目标基因相关的互作数据
     */
    @Select("SELECT * FROM tf WHERE target_gene = #{targetGene}")
    List<TFData> selectTFDataByTargetGene(@Param("targetGene") String targetGene);
    
    /**
     * 查询与特定基因相关的互作数据（作为转录因子或目标基因）
     */
    @Select("SELECT * FROM tf WHERE tf = #{gene} OR target_gene = #{gene}")
    List<TFData> selectTFDataByGene(@Param("gene") String gene);
}

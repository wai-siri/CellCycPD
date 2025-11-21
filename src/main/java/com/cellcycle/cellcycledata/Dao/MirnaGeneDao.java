package com.cellcycle.cellcycledata.Dao;

import org.apache.ibatis.annotations.Mapper;
import org.apache.ibatis.annotations.Param;
import org.apache.ibatis.annotations.Select;
import com.cellcycle.cellcycledata.Entity.MirnaGene;

import java.util.List;
import java.util.Map;

/**
 * MirnaGene数据访问接口
 */
@Mapper
public interface MirnaGeneDao {
    /**
     * 查询所有miRNA-gene互作数据
     */
    @Select("SELECT * FROM mirna_gene")
    List<MirnaGene> selectAllMirnaGeneData();
    
    /**
     * 调试：查询表结构
     */
    @Select("SHOW COLUMNS FROM mirna_gene")
    List<Map<String, Object>> showColumns();
    
    /**
     * 调试：原始SQL查询
     */
    @Select("SELECT miRNA, gene FROM mirna_gene LIMIT 1")
    List<Map<String, Object>> selectRawData();
    
    /**
     * 根据miRNA名称查询互作数据
     */
    @Select("SELECT * FROM mirna_gene WHERE miRNA = #{miRNA}")
    List<MirnaGene> selectMirnaGeneDataByMirna(@Param("miRNA") String miRNA);
    
    /**
     * 根据基因名称查询互作数据
     */
    @Select("SELECT * FROM mirna_gene WHERE gene = #{gene}")
    List<MirnaGene> selectMirnaGeneDataByGene(@Param("gene") String gene);
    
    /**
     * 查询所有唯一的miRNA名称
     */
    @Select("SELECT DISTINCT miRNA FROM mirna_gene ORDER BY miRNA")
    List<String> selectAllDistinctMirnas();
    
    /**
     * 查询所有唯一的基因名称
     */
    @Select("SELECT DISTINCT gene FROM mirna_gene ORDER BY gene")
    List<String> selectAllDistinctGenes();
    
    /**
     * 根据miRNA或基因名称查询互作数据
     */
    @Select("SELECT * FROM mirna_gene WHERE miRNA = #{name} OR gene = #{name}")
    List<MirnaGene> selectMirnaGeneDataByMirnaOrGene(@Param("name") String name);
}

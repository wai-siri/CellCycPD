package com.cellcycle.cellcycledata.Service;

import org.springframework.stereotype.Service;
import javax.annotation.Resource;
import com.cellcycle.cellcycledata.Dao.MirnaGeneDao;
import com.cellcycle.cellcycledata.Entity.MirnaGene;

import java.util.List;
import java.util.Map;

/**
 * MirnaGene服务类，处理miRNA-gene互作数据的业务逻辑
 */
@Service("mirnaGeneService")
public class MirnaGeneService {
    @Resource
    private MirnaGeneDao mirnaGeneDao;
    
    /**
     * 获取所有miRNA-gene互作数据
     */
    public List<MirnaGene> getAllMirnaGeneData() {
        // 调试：查看表结构
        try {
            System.out.println("=== Table Structure Debug ===");
            List<Map<String, Object>> columns = this.mirnaGeneDao.showColumns();
            for (Map<String, Object> column : columns) {
                System.out.println("Column: " + column);
            }
            
            System.out.println("=== Raw Data Debug ===");
            List<Map<String, Object>> rawData = this.mirnaGeneDao.selectRawData();
            for (Map<String, Object> row : rawData) {
                System.out.println("Raw row: " + row);
            }
            System.out.println("=== End Raw Debug ===");
        } catch (Exception e) {
            System.out.println("Debug query failed: " + e.getMessage());
        }
        
        return this.mirnaGeneDao.selectAllMirnaGeneData();
    }
    
    /**
     * 根据miRNA名称获取互作数据
     */
    public List<MirnaGene> getMirnaGeneDataByMirna(String miRNA) {
        return this.mirnaGeneDao.selectMirnaGeneDataByMirna(miRNA);
    }
    
    /**
     * 根据基因名称获取互作数据
     */
    public List<MirnaGene> getMirnaGeneDataByGene(String gene) {
        return this.mirnaGeneDao.selectMirnaGeneDataByGene(gene);
    }
    
    /**
     * 获取所有唯一的miRNA名称
     */
    public List<String> getAllDistinctMirnas() {
        return this.mirnaGeneDao.selectAllDistinctMirnas();
    }
    
    /**
     * 获取所有唯一的基因名称
     */
    public List<String> getAllDistinctGenes() {
        return this.mirnaGeneDao.selectAllDistinctGenes();
    }
    
    /**
     * 根据miRNA或基因名称获取互作数据
     */
    public List<MirnaGene> getMirnaGeneDataByMirnaOrGene(String name) {
        return this.mirnaGeneDao.selectMirnaGeneDataByMirnaOrGene(name);
    }
}

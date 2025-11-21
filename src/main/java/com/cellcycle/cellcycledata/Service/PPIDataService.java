package com.cellcycle.cellcycledata.Service;

import org.springframework.stereotype.Service;
import javax.annotation.Resource;
import com.cellcycle.cellcycledata.Dao.PPIDataDao;
import com.cellcycle.cellcycledata.Entity.PPIData;

import java.util.List;
import java.util.Map;

@Service("ppiDataService")
public class PPIDataService {
    @Resource
    private PPIDataDao ppiDataDao;
    
    /**
     * 获取所有蛋白质互作数据
     */
    public List<PPIData> getAllPPIData() {
        // 调试：查看表结构
        try {
            System.out.println("=== Table Structure Debug ===");
            List<Map<String, Object>> columns = this.ppiDataDao.showColumns();
            for (Map<String, Object> column : columns) {
                System.out.println("Column: " + column);
            }
            
            System.out.println("=== Raw Data Debug ===");
            List<Map<String, Object>> rawData = this.ppiDataDao.selectRawData();
            for (Map<String, Object> row : rawData) {
                System.out.println("Raw row: " + row);
            }
            System.out.println("=== End Raw Debug ===");
        } catch (Exception e) {
            System.out.println("Debug query failed: " + e.getMessage());
        }
        
        return this.ppiDataDao.selectAllPPIData();
    }
    
    /**
     * 根据最低分数阈值获取蛋白质互作数据
     */
    public List<PPIData> getPPIDataByMinScore(Double minScore) {
        return this.ppiDataDao.selectPPIDataByMinScore(minScore);
    }
    
    /**
     * 获取与特定基因相关的蛋白质互作数据
     */
    public List<PPIData> getPPIDataByGene(String gene) {
        return this.ppiDataDao.selectPPIDataByGene(gene);
    }
    
    /**
     * 根据基因和最低分数阈值获取蛋白质互作数据
     */
    public List<PPIData> getPPIDataByGeneAndMinScore(String gene, Double minScore) {
        return this.ppiDataDao.selectPPIDataByGeneAndMinScore(gene, minScore);
    }
}

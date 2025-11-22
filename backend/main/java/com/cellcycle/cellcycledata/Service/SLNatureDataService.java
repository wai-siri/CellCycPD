package com.cellcycle.cellcycledata.Service;

import org.springframework.stereotype.Service;
import javax.annotation.Resource;
import com.cellcycle.cellcycledata.Dao.SLNatureDataDao;
import com.cellcycle.cellcycledata.Entity.SLNatureData;

import java.util.List;
import java.util.Map;

@Service("slNatureDataService")
public class SLNatureDataService {
    @Resource
    private SLNatureDataDao slNatureDataDao;
    
    /**
     * 获取所有合成致死自然数据
     */
    public List<SLNatureData> selectAllSLNatureData() {
        // 调试：查看表结构
        try {
            System.out.println("=== SL Nature Table Structure Debug ===");
            List<Map<String, Object>> columns = this.slNatureDataDao.showColumns();
            for (Map<String, Object> column : columns) {
                System.out.println("Column: " + column);
            }
            
            System.out.println("=== SL Nature Raw Data Debug ===");
            List<Map<String, Object>> rawData = this.slNatureDataDao.selectRawData();
            for (Map<String, Object> row : rawData) {
                System.out.println("Raw row: " + row);
            }
            System.out.println("=== End SL Nature Raw Debug ===");
        } catch (Exception e) {
            System.out.println("Debug query failed: " + e.getMessage());
        }
        
        // 获取所有数据并打印详细日志
        List<SLNatureData> allData = this.slNatureDataDao.selectAllSLNatureData();
        System.out.println("=== SL Nature All Data Debug ===");
        System.out.println("Total records: " + allData.size());
        if (!allData.isEmpty()) {
            SLNatureData firstRecord = allData.get(0);
            System.out.println("First record details:");
            System.out.println("  geneA: " + firstRecord.getGeneA());
            System.out.println("  stageGeneA: " + firstRecord.getStageGeneA());
            System.out.println("  geneB: " + firstRecord.getGeneB());
            System.out.println("  stageGeneB: " + firstRecord.getStageGeneB());
            System.out.println("  geminiSensitive: " + firstRecord.getGeminiSensitive());
            System.out.println("  cellLine: " + firstRecord.getCellLine());
        }
        System.out.println("=== End SL Nature All Data Debug ===");
        
        return allData;
    }
    
    /**
     * 根据基因名称查询相关的合成致死自然数据
     */
    public List<SLNatureData> selectSLNatureDataByGene(String gene) {
        return this.slNatureDataDao.selectSLNatureDataByGene(gene);
    }
}

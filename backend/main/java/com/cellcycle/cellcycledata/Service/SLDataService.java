package com.cellcycle.cellcycledata.Service;

import org.springframework.stereotype.Service;
import javax.annotation.Resource;
import com.cellcycle.cellcycledata.Dao.SLDataDao;
import com.cellcycle.cellcycledata.Entity.SLData;

import java.util.List;
import java.util.Map;

@Service("slDataService")
public class SLDataService {
    @Resource
    private SLDataDao slDataDao;
    
    /**
     * 获取所有合成致死数据
     */
    public List<SLData> selectAllSLData() {
        // 调试：查看表结构
        try {
            System.out.println("=== SL Table Structure Debug ===");
            List<Map<String, Object>> columns = this.slDataDao.showColumns();
            for (Map<String, Object> column : columns) {
                System.out.println("Column: " + column);
            }
            
            System.out.println("=== SL Raw Data Debug ===");
            List<Map<String, Object>> rawData = this.slDataDao.selectRawData();
            for (Map<String, Object> row : rawData) {
                System.out.println("Raw row: " + row);
            }
            System.out.println("=== End SL Raw Debug ===");
        } catch (Exception e) {
            System.out.println("Debug query failed: " + e.getMessage());
        }
        
        return this.slDataDao.selectAllSLData();
    }
    
    /**
     * 根据基因名称查询相关的合成致死数据
     */
    public List<SLData> selectSLDataByGene(String gene) {
        return this.slDataDao.selectSLDataByGene(gene);
    }
}

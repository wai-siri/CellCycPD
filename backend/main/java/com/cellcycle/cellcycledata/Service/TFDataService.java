package com.cellcycle.cellcycledata.Service;

import org.springframework.stereotype.Service;
import javax.annotation.Resource;
import com.cellcycle.cellcycledata.Dao.TFDataDao;
import com.cellcycle.cellcycledata.Entity.TFData;

import java.util.List;
import java.util.Map;

/**
 * 转录因子-基因互作数据服务类
 */
@Service("tfDataService")
public class TFDataService {
    @Resource
    private TFDataDao tfDataDao;
    
    /**
     * 获取所有转录因子-基因互作数据
     */
    public List<TFData> getAllTFData() {
        // 调试：查看表结构
        try {
            System.out.println("=== TF Table Structure Debug ===");
            List<Map<String, Object>> columns = this.tfDataDao.showColumns();
            for (Map<String, Object> column : columns) {
                System.out.println("Column: " + column);
            }
            
            System.out.println("=== TF Raw Data Debug ===");
            List<Map<String, Object>> rawData = this.tfDataDao.selectRawData();
            for (Map<String, Object> row : rawData) {
                System.out.println("Raw row: " + row);
            }
            System.out.println("=== End TF Raw Debug ===");
        } catch (Exception e) {
            System.out.println("TF Debug query failed: " + e.getMessage());
        }
        
        return this.tfDataDao.selectAllTFData();
    }
    
    /**
     * 获取与特定转录因子相关的互作数据
     */
    public List<TFData> getTFDataByTF(String tf) {
        return this.tfDataDao.selectTFDataByTF(tf);
    }
    
    /**
     * 获取与特定目标基因相关的互作数据
     */
    public List<TFData> getTFDataByTargetGene(String targetGene) {
        return this.tfDataDao.selectTFDataByTargetGene(targetGene);
    }
    
    /**
     * 获取与特定基因相关的互作数据（作为转录因子或目标基因）
     */
    public List<TFData> getTFDataByGene(String gene) {
        return this.tfDataDao.selectTFDataByGene(gene);
    }
}

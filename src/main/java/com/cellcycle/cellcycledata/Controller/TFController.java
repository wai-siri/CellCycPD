package com.cellcycle.cellcycledata.Controller;

import com.cellcycle.cellcycledata.Entity.TFData;
import com.cellcycle.cellcycledata.Service.TFDataService;
import org.springframework.web.bind.annotation.*;

import javax.annotation.Resource;
import java.util.List;

/**
 * 转录因子-基因互作数据控制器
 */
@RestController
@CrossOrigin
@RequestMapping("cellcycledata")
@ResponseBody
public class TFController {
    
    @Resource
    private TFDataService tfDataService;
    
    /**
     * 获取所有转录因子-基因互作数据
     */
    @GetMapping("/tf/all")
    public List<TFData> getAllTFData() {
        List<TFData> result = this.tfDataService.getAllTFData();
        
        // 调试：打印第一条数据的所有字段
        if (!result.isEmpty()) {
            TFData first = result.get(0);
            System.out.println("=== Backend Debug: First TF Data ===");
            System.out.println("tf: " + first.getTf());
            System.out.println("target_gene: " + first.getTarget_gene());
            System.out.println("interaction: " + first.getInteraction());
            System.out.println("=== End Backend Debug ===");
        }
        
        return result;
    }
    
    /**
     * 获取与特定转录因子相关的互作数据
     */
    @GetMapping("/tf/factor/{tf}")
    public List<TFData> getTFDataByTF(@PathVariable("tf") String tf) {
        return this.tfDataService.getTFDataByTF(tf);
    }
    
    /**
     * 获取与特定目标基因相关的互作数据
     */
    @GetMapping("/tf/target/{targetGene}")
    public List<TFData> getTFDataByTargetGene(@PathVariable("targetGene") String targetGene) {
        return this.tfDataService.getTFDataByTargetGene(targetGene);
    }
    
    /**
     * 获取与特定基因相关的互作数据（作为转录因子或目标基因）
     */
    @GetMapping("/tf/gene/{gene}")
    public List<TFData> getTFDataByGene(@PathVariable("gene") String gene) {
        return this.tfDataService.getTFDataByGene(gene);
    }
}

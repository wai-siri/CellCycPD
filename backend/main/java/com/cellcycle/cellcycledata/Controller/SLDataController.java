package com.cellcycle.cellcycledata.Controller;

import com.cellcycle.cellcycledata.Entity.SLData;
import com.cellcycle.cellcycledata.Service.SLDataService;
import org.springframework.web.bind.annotation.*;

import javax.annotation.Resource;
import java.util.List;

/**
 * 合成致死(Synthetic Lethality)数据控制器
 */
@RestController
@CrossOrigin
@RequestMapping("cellcycledata")
@ResponseBody
public class SLDataController {
    
    @Resource
    private SLDataService slDataService;
    
    /**
     * 获取所有合成致死数据
     * @return 所有合成致死数据列表
     */
    @GetMapping("/sl/all")
    public List<SLData> getAllSLData() {
        try {
            // 直接调用服务层方法，其中已包含调试代码
            List<SLData> result = slDataService.selectAllSLData();
            
            // 返回结果
            return result;
        } catch (Exception e) {
            System.err.println("Error in getAllSLData: " + e.getMessage());
            e.printStackTrace();
            throw e;
        }
    }
    
    /**
     * 根据基因名称获取相关的合成致死数据
     * @param gene 基因名称
     * @return 与指定基因相关的合成致死数据列表
     */
    @GetMapping("/sl/gene/{gene}")
    public List<SLData> getSLDataByGene(@PathVariable("gene") String gene) {
        return slDataService.selectSLDataByGene(gene);
    }
}

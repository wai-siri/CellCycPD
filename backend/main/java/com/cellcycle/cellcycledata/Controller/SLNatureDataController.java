package com.cellcycle.cellcycledata.Controller;

import com.cellcycle.cellcycledata.Entity.SLNatureData;
import com.cellcycle.cellcycledata.Service.SLNatureDataService;
import org.springframework.web.bind.annotation.*;

import javax.annotation.Resource;
import java.util.List;

/**
 * 合成致死自然数据(SL Nature)控制器
 */
@RestController
@CrossOrigin
@RequestMapping("cellcycledata")
@ResponseBody
public class SLNatureDataController {
    
    @Resource
    private SLNatureDataService slNatureDataService;
    
    /**
     * 获取所有合成致死自然数据
     * @return 所有合成致死自然数据列表
     */
    @GetMapping("/sl/nature/all")
    public List<SLNatureData> getAllSLNatureData() {
        try {
            // 直接调用服务层方法，其中已包含调试代码
            List<SLNatureData> result = slNatureDataService.selectAllSLNatureData();
            
            // 返回结果
            return result;
        } catch (Exception e) {
            System.err.println("Error in getAllSLNatureData: " + e.getMessage());
            e.printStackTrace();
            throw e;
        }
    }
    
    /**
     * 根据基因名称获取相关的合成致死自然数据
     * @param gene 基因名称
     * @return 与指定基因相关的合成致死自然数据列表
     */
    @GetMapping("/sl/nature/gene/{gene}")
    public List<SLNatureData> getSLNatureDataByGene(@PathVariable("gene") String gene) {
        return slNatureDataService.selectSLNatureDataByGene(gene);
    }
}

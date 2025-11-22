package com.cellcycle.cellcycledata.Controller;

import com.cellcycle.cellcycledata.Entity.PPIData;
import com.cellcycle.cellcycledata.Service.PPIDataService;
import org.springframework.web.bind.annotation.*;

import javax.annotation.Resource;
import java.util.List;

@RestController
@CrossOrigin
@RequestMapping("cellcycledata")
@ResponseBody
public class PPIController {
    
    @Resource
    private PPIDataService ppiDataService;
    
    /**
     * 获取所有蛋白质互作数据
     */
    @GetMapping("/ppi/all")
    public List<PPIData> getAllPPIData() {
        List<PPIData> result = this.ppiDataService.getAllPPIData();
        
        // 调试：打印第一条数据的所有字段
        if (!result.isEmpty()) {
            PPIData first = result.get(0);
            System.out.println("=== Backend Debug: First PPI Data ===");
            System.out.println("node1: " + first.getNode1());
            System.out.println("node2: " + first.getNode2());
            System.out.println("node1_string_id: " + first.getNode1_string_id());
            System.out.println("node2_string_id: " + first.getNode2_string_id());
            System.out.println("neighborhood_on_chromosome: " + first.getNeighborhood_on_chromosome());
            System.out.println("gene_fusion: " + first.getGene_fusion());
            System.out.println("phylogenetic_cooccurrence: " + first.getPhylogenetic_cooccurrence());
            System.out.println("homology: " + first.getHomology());
            System.out.println("coexpression: " + first.getCoexpression());
            System.out.println("experimentally_determined_interaction: " + first.getExperimentally_determined_interaction());
            System.out.println("database_annotated: " + first.getDatabase_annotated());
            System.out.println("automated_textmining: " + first.getAutomated_textmining());
            System.out.println("combined_score: " + first.getCombined_score());
            System.out.println("=== End Backend Debug ===");
        }
        
        return result;
    }
    
    /**
     * 根据最低分数阈值获取蛋白质互作数据
     */
    @GetMapping("/ppi/score/{minScore}")
    public List<PPIData> getPPIDataByMinScore(@PathVariable("minScore") Double minScore) {
        return this.ppiDataService.getPPIDataByMinScore(minScore);
    }
    
    /**
     * 获取与特定基因相关的蛋白质互作数据
     */
    @GetMapping("/ppi/gene/{gene}")
    public List<PPIData> getPPIDataByGene(@PathVariable("gene") String gene) {
        return this.ppiDataService.getPPIDataByGene(gene);
    }
    
    /**
     * 根据基因和最低分数阈值获取蛋白质互作数据
     */
    @GetMapping("/ppi/gene/{gene}/score/{minScore}")
    public List<PPIData> getPPIDataByGeneAndMinScore(
            @PathVariable("gene") String gene,
            @PathVariable("minScore") Double minScore) {
        return this.ppiDataService.getPPIDataByGeneAndMinScore(gene, minScore);
    }
}

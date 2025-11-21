package com.cellcycle.cellcycledata.Controller;

import com.cellcycle.cellcycledata.Entity.GeneTarget;
import com.cellcycle.cellcycledata.Service.GeneTargetService;
import javax.annotation.Resource;
import org.springframework.web.bind.annotation.*;

import java.util.List;

@RestController
@CrossOrigin
@RequestMapping("cellcycledata")
@ResponseBody
public class GeneTargetController {
    
    @Resource
    private GeneTargetService geneTargetService;

    @GetMapping("/gene-targets")
    public List<GeneTarget> getAllGeneTargets() {
        return this.geneTargetService.selectAllGeneTargets();
    }

    @GetMapping("/gene-targets/gene/{gene}")
    public List<GeneTarget> getGeneTargetsByGene(@PathVariable String gene) {
        return this.geneTargetService.selectByGene(gene);
    }

    @GetMapping("/gene-targets/stage/{stage}")
    public List<GeneTarget> getGeneTargetsByStage(@PathVariable String stage) {
        return this.geneTargetService.selectByStage(stage);
    }
}

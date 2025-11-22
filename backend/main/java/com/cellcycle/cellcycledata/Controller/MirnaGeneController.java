package com.cellcycle.cellcycledata.Controller;

import com.cellcycle.cellcycledata.Entity.MirnaGene;
import com.cellcycle.cellcycledata.Service.MirnaGeneService;
import org.springframework.web.bind.annotation.*;

import javax.annotation.Resource;
import java.util.List;

/**
 * MirnaGene控制器，提供miRNA-gene互作数据的REST API接口
 */
@RestController
@RequestMapping("cellcycledata")
public class MirnaGeneController {

    @Resource
    private MirnaGeneService mirnaGeneService;

    /**
     * 获取所有miRNA-gene互作数据
     */
    @GetMapping("/mirna/all")
    public List<MirnaGene> getAllMirnaGeneData() {
        List<MirnaGene> result = this.mirnaGeneService.getAllMirnaGeneData();
        
        // 调试：打印第一条数据的所有字段
        if (!result.isEmpty()) {
            MirnaGene first = result.get(0);
            System.out.println("=== Backend Debug: First MirnaGene Data ===");
            System.out.println("miRNA: " + first.getMiRNA());
            System.out.println("gene: " + first.getGene());
            System.out.println("=== End Backend Debug ===");
        }
        
        return result;
    }

    /**
     * 根据miRNA名称获取互作数据
     */
    @GetMapping("/mirna/mirna/{miRNA}")
    public List<MirnaGene> getMirnaGeneDataByMirna(@PathVariable("miRNA") String miRNA) {
        return this.mirnaGeneService.getMirnaGeneDataByMirna(miRNA);
    }

    /**
     * 根据基因名称获取互作数据
     */
    @GetMapping("/mirna/gene/{gene}")
    public List<MirnaGene> getMirnaGeneDataByGene(@PathVariable("gene") String gene) {
        return this.mirnaGeneService.getMirnaGeneDataByGene(gene);
    }

    /**
     * 获取所有唯一的miRNA名称
     */
    @GetMapping("/mirna/mirnas")
    public List<String> getAllDistinctMirnas() {
        return this.mirnaGeneService.getAllDistinctMirnas();
    }

    /**
     * 获取所有唯一的基因名称
     */
    @GetMapping("/mirna/genes")
    public List<String> getAllDistinctGenes() {
        return this.mirnaGeneService.getAllDistinctGenes();
    }

    /**
     * 根据miRNA或基因名称获取互作数据
     */
    @GetMapping("/mirna/search/{name}")
    public List<MirnaGene> getMirnaGeneDataByMirnaOrGene(@PathVariable("name") String name) {
        return this.mirnaGeneService.getMirnaGeneDataByMirnaOrGene(name);
    }
}

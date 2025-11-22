package com.cellcycle.cellcycledata.Controller;

import com.cellcycle.cellcycledata.Entity.cellcycledata;
import com.cellcycle.cellcycledata.Service.cellcycledataService;
import javax.annotation.Resource;
import org.springframework.web.bind.annotation.*;
import org.springframework.core.io.ClassPathResource;
import org.springframework.http.ResponseEntity;
import org.springframework.http.HttpStatus;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

@RestController
@CrossOrigin
@RequestMapping("cellcycledata")
@ResponseBody
public class cellcycledataController {
    
    @Resource
    private cellcycledataService cellcycledataService;

    @GetMapping("/allcellcyclegene")
    public List<cellcycledata> selectAllCellcycledata() {
        return this.cellcycledataService.selectAllCellcycledata();
    }

    @GetMapping("/allcellcyclegenename")
    public String[] GetAllGeneName(){
        List<cellcycledata> cellcycledata = this.cellcycledataService.selectAllCellcycledata();
        String[] gene = new String[cellcycledata.size()];
        for (int i = 0; i < cellcycledata.size(); i++) {
            gene[i] = cellcycledata.get(i).getGene_Symbol();
        }
        return gene;
    }

    @GetMapping("/kegg-pathway")
    public ResponseEntity<?> getKEGGPathway() {
        try {
            ClassPathResource resource = new ClassPathResource("static/Browser/KEGG_pathway.txt");
            List<Map<String, String>> results = new ArrayList<>();
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(resource.getInputStream()))) {
                String line;
                String[] headers = null;
                while ((line = reader.readLine()) != null) {
                    String[] values = line.split("\t");
                    if (headers == null) {
                        headers = values;
                    } else {
                        Map<String, String> row = new HashMap<>();
                        for (int i = 0; i < headers.length && i < values.length; i++) {
                            row.put(headers[i], values[i]);
                        }
                        results.add(row);
                    }
                }
            }
            return ResponseEntity.ok(results);
        } catch (IOException e) {
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Error reading KEGG pathway data");
        }
    }
    
    @GetMapping("/hallmark-cell")
    public ResponseEntity<?> getHallmarkCellData() {
        try {
            ClassPathResource resource = new ClassPathResource("static/home/hallmark_cell.txt");
            List<Map<String, Object>> results = new ArrayList<>();
            Map<String, List<String>> geneHallmarkMap = new HashMap<>();
            Map<String, Integer> geneCountMap = new HashMap<>();
            
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(resource.getInputStream()))) {
                String line;
                String[] headers = null;
                while ((line = reader.readLine()) != null) {
                    String[] values = line.split("\t");
                    if (headers == null) {
                        headers = values;
                    } else if (values.length >= 2) {
                        String gene = values[0];
                        String hallmark = values[1];
                        
                        // 统计基因出现次数
                        geneCountMap.put(gene, geneCountMap.getOrDefault(gene, 0) + 1);
                        
                        // 收集基因对应的所有hallmark
                        if (!geneHallmarkMap.containsKey(gene)) {
                            geneHallmarkMap.put(gene, new ArrayList<>());
                        }
                        geneHallmarkMap.get(gene).add(hallmark);
                    }
                }
            }
            
            // 构建结果
            for (String gene : geneHallmarkMap.keySet()) {
                Map<String, Object> row = new HashMap<>();
                row.put("gene", gene);
                row.put("hallmarks", geneHallmarkMap.get(gene));
                row.put("count", geneCountMap.get(gene));
                results.add(row);
            }
            
            return ResponseEntity.ok(results);
        } catch (IOException e) {
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Error reading hallmark cell data");
        }
    }
    
    @GetMapping("/mirna-gene-wordcloud")
    public ResponseEntity<?> getMirnaGeneWordcloud() {
        try {
            ClassPathResource resource = new ClassPathResource("static/network/miRNA_gene_interaction_number.txt");
            List<Map<String, Object>> results = new ArrayList<>();
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(resource.getInputStream()))) {
                String line;
                String[] headers = null;
                while ((line = reader.readLine()) != null) {
                    String[] values = line.split("\t");
                    if (headers == null) {
                        headers = values;
                    } else if (values.length >= 2) {
                        Map<String, Object> row = new HashMap<>();
                        row.put("name", values[0]); // 基因名称
                        try {
                            row.put("value", Integer.parseInt(values[1])); // 频数
                        } catch (NumberFormatException e) {
                            row.put("value", 0);
                        }
                        results.add(row);
                    }
                }
            }
            return ResponseEntity.ok(results);
        } catch (IOException e) {
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Error reading miRNA-gene wordcloud data");
        }
    }
} 
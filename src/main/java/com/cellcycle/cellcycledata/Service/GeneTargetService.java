package com.cellcycle.cellcycledata.Service;

import org.springframework.stereotype.Service;
import javax.annotation.Resource;
import com.cellcycle.cellcycledata.Dao.GeneTargetDao;
import com.cellcycle.cellcycledata.Entity.GeneTarget;

import java.util.List;

@Service("geneTargetService")
public class GeneTargetService {
    @Resource
    private GeneTargetDao geneTargetDao;
    
    public List<GeneTarget> selectAllGeneTargets() {
        return this.geneTargetDao.selectAllGeneTargets();
    }
    
    public List<GeneTarget> selectByGene(String gene) {
        return this.geneTargetDao.selectByGene(gene);
    }
    
    public List<GeneTarget> selectByStage(String stage) {
        return this.geneTargetDao.selectByStage(stage);
    }
}

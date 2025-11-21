package com.cellcycle.cellcycledata.Controller;

import com.cellcycle.cellcycledata.Entity.DrugInformation;
import com.cellcycle.cellcycledata.Service.DrugInformationService;
import javax.annotation.Resource;
import org.springframework.web.bind.annotation.*;

import java.util.List;

@RestController
@CrossOrigin
@RequestMapping("cellcycledata")
@ResponseBody
public class DrugInformationController {
    
    @Resource
    private DrugInformationService drugInformationService;

    @GetMapping("/drug-information")
    public List<DrugInformation> getAllDrugInformation() {
        return this.drugInformationService.selectAllDrugInformation();
    }

    @GetMapping("/drug-information/{drug}")
    public DrugInformation getDrugInformationByDrug(@PathVariable String drug) {
        return this.drugInformationService.selectByDrug(drug);
    }
}

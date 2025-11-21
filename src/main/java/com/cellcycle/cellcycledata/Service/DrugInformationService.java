package com.cellcycle.cellcycledata.Service;

import org.springframework.stereotype.Service;
import javax.annotation.Resource;
import com.cellcycle.cellcycledata.Dao.DrugInformationDao;
import com.cellcycle.cellcycledata.Entity.DrugInformation;

import java.util.List;

@Service("drugInformationService")
public class DrugInformationService {
    @Resource
    private DrugInformationDao drugInformationDao;
    
    public List<DrugInformation> selectAllDrugInformation() {
        return this.drugInformationDao.selectAllDrugInformation();
    }
    
    public DrugInformation selectByDrug(String drug) {
        return this.drugInformationDao.selectByDrug(drug);
    }
}

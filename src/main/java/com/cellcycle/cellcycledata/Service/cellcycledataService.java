package com.cellcycle.cellcycledata.Service;

import org.springframework.stereotype.Service;
import javax.annotation.Resource;
import com.cellcycle.cellcycledata.Dao.cellcycledataDao;
import com.cellcycle.cellcycledata.Entity.cellcycledata;

import java.util.List;

@Service("cellcycledataService")

public class cellcycledataService  {
    @Resource
    private cellcycledataDao cellcycledataDao;
  
  
    public List<cellcycledata> selectAllCellcycledata(){return this.cellcycledataDao.selectAllCellcycledata();}
  
  
  
  }
package com.cellcycle.cellcycledata.Dao;

import org.apache.ibatis.annotations.Mapper;
import com.cellcycle.cellcycledata.Entity.cellcycledata;

import java.util.List;

@Mapper
public interface cellcycledataDao {

    //  查询表格中所有睡眠相关基因
      List<cellcycledata>selectAllCellcycledata();
    
    
    }
    
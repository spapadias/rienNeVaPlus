-- Script for accuracy analysis 

-- Delete entries from all tables
delete from ls_sc0_speed;
-- delete from ls_sc1_speed;
-- delete from ls_sc2_speed;
delete from ls_sc3_speed;
-- delete from ls_sc4_speed;
-- delete from ls_sc5_speed;	

-- Real - Rome
COPY LS_sc0_speed FROM 'C:\Program Files\PostgreSQL\s0.txt'  DELIMITER ',' NULL ''; -- NO_LS
-- COPY LS_sc1_speed FROM 'C:\Program Files\PostgreSQL\s1.txt'  DELIMITER ',' NULL ''; -- RND
-- COPY LS_sc2_speed FROM 'C:\Program Files\PostgreSQL\s2.txt'  DELIMITER ',' NULL ''; -- INV
COPY LS_sc3_speed FROM 'C:\Program Files\PostgreSQL\s3.txt'  DELIMITER ',' NULL ''; -- COMB
-- COPY LS_sc4_speed FROM 'C:\Program Files\PostgreSQL\s4.txt'  DELIMITER ',' NULL ''; -- 
-- COPY LS_sc5_speed FROM 'C:\Program Files\PostgreSQL\s5.txt'  DELIMITER ',' NULL ''; -- 

 
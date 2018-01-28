-- Script for accuracy analysis 
 /*
-- Delete entries from all tables
delete from ls_sc0_speed;
delete from ls_sc1_speed;
delete from ls_sc2_speed;

-- Real - Rome
COPY LS_sc0_speed FROM 'C:\Program Files\PostgreSQL\s0.txt'  DELIMITER ',' NULL ''; -- NO_LS
COPY LS_sc1_speed FROM 'C:\Program Files\PostgreSQL\s1.txt'  DELIMITER ',' NULL ''; -- R
COPY LS_sc2_speed FROM 'C:\Program Files\PostgreSQL\s2.txt'  DELIMITER ',' NULL ''; -- TF_now (no history)
 */ 
-- /*
-- Mean and standard deviation among all difference values in speed estimates
WITH speed_differences AS 
(
    SELECT 
        A.t, 
        A.qid,
        A.speed                AS est_Speed_LS, 
        B.speed                AS real_Speed_noLS, 
        abs(A.speed - B.speed) AS diff_Speed,
    	-- (abs(A.speed - B.speed) / B.speed) AS relative_diff_Speed, 
        A.num_objects          AS used_Objects_LS, 
        B.num_objects          AS all_Objects_noLS

    FROM
    -- *********************************** --
        LS_sc0_speed A, LS_sc0_speed B
    -- *********************************** --

    WHERE 
        A.t      = B.t   AND 
        A.qID    = B.qID AND 
        A.speed >= 0     AND 
        B.speed >= 0
        -- CAUTION: Include only the execution cycles where LS took place

    	-- For capacity=10% its from e.c. t=1240 till t=2970 
    	-- AND A.t >= 1500 AND B.t <= 3540
)
SELECT 
        round(AVG(diff_speed)::numeric,3)    AS avg_Diff, 
        round(stddev(diff_speed)::numeric,3) AS stdev_avg_Diff 
        -- round(AVG(relative_diff_speed)::numeric,3) AS rel_avg_Diff,
        -- round(stddev(relative_diff_speed)::numeric,3) AS stdev_rel_avg_Diff
FROM speed_differences;
-- */
-- End of script
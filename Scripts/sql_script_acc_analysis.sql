-- Script for accuracy analysis 

-- Delete entries from all tables
delete from ls_sc0_speed;
delete from ls_sc1_speed;
delete from ls_sc2_speed;
delete from ls_sc3_speed;

-- Insert entries to all tables
COPY LS_sc0_speed FROM 'C:\sc0.txt'  DELIMITER ',' NULL ''; -- NO_LS
COPY LS_sc1_speed FROM 'C:\sc1.txt'  DELIMITER ',' NULL ''; -- R
COPY LS_sc2_speed FROM 'C:\sc2.txt'  DELIMITER ',' NULL ''; -- TF
COPY LS_sc3_speed FROM 'C:\sc3.txt'  DELIMITER ',' NULL ''; -- TFw

-- Mean and standard deviation among all difference values in speed estimates
WITH speed_differences AS 
(
    SELECT 
        A.t, 
        A.qid,
        A.speed                AS est_Speed_LS, 
        B.speed                AS real_Speed_noLS, 
        abs(A.speed - B.speed) AS diff_Speed, -- (abs(A.speed - B.speed) / B.speed) AS relative_diff_Speed, 
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
        -- For capacity=50% 
        -- AND A.t >= 1240 AND B.t <= 2970
    	-- For capacity=10% its from e.c. t=1240 till t=2970 
    	-- AND A.t >= 850 AND B.t <= 3360
)
SELECT 
        round(AVG(diff_speed)::numeric,3)    AS avg_Diff, -- round(AVG(relative_diff_speed)::numeric,3) AS rel_avg_Diff,
        round(stddev(diff_speed)::numeric,3) AS stdev_avg_Diff -- round(stddev(relative_diff_speed)::numeric,3) AS stdev_rel_avg_Diff
FROM speed_differences;

-- End of script



-- Big speed differences
/*
SELECT 
	A.t, A.qid, 
    A.speed AS est_Speed_LS, 
    B.speed AS real_Speed_noLS, 
    round(abs(A.speed - B.speed)::numeric,3) AS diff_Speed, 
    round((abs(A.speed - B.speed)/ B.speed)::numeric,3)  AS rel_diff_Speed, 
    A.num_objects AS used_Objects_LS, 
    B.num_objects AS all_Objects_noLS
FROM
	-- Insert the experiment for comparison with ground truth here
	LS_sc1_speed A, LS_sc0_speed B
    -- ***********************************************************
WHERE 
	A.t = B.t
	AND A.qID = B.qID
	AND A.speed >0 AND B.speed>0
	AND abs(A.speed - B.speed)>20
	-- CAUTION: Include only the execution cycles where LS took place
    -- For capacity=50% its from e.c. t=1240 till t=2970 
    -- AND A.t >= 1240 AND B.t <= 2970
    -- For capacity=50% its from e.c. t=1240 till t=2970 
    AND A.t >= 850 AND B.t <= 3360
ORDER BY 
	t ASC, abs(A.speed - B.speed) DESC;
*/
-- Speed differences w.r.t. axes type
/*
WITH 
speed_differences AS 
(
	SELECT 
    	A.t, A.qid, 
        A.speed AS est_Speed_LS, 
    	B.speed AS real_Speed_noLS, 
    	abs(A.speed - B.speed) AS diff_Speed, 
    	A.num_objects AS used_Objects_LS, 
    	B.num_objects AS all_Objects_noLS
    FROM
    	-- ************************************************* --
    	LS_sc1_speed A, LS_sc0_speed B
	    -- ************************************************* --
    WHERE 
    	A.t = B.t
        -- AND A.t >=1200
    	AND A.qID = B.qID
    	AND A.speed >0 AND B.speed>0
        -- CAUTION: Include only the execution cycles where LS took place
        -- For capacity=50% its from e.c. t=1240 till t=2970 
        AND A.t >= 1240 
        AND B.t <= 2970
    	-- CAUTION: Include only the execution cycles where LS took place
        -- For capacity=50% its from e.c. t=1240 till t=2970 
        AND A.t >= 1240 
        AND B.t <= 2970
),
axes_classification AS 
(
    SELECT 
    	qid,
    	min(win_slide) AS class
    FROM 
    	-- Insert the query file here
    	athens_axes
    	-- Query file inserted
    GROUP BY 
    	qid
),
queries_per_class AS 
(
    SELECT 
    	win_range, 
    	win_slide, 
    	count (DISTINCT qid) AS num_queries
    FROM 
    	-- Insert the query file here
    	athens_axes
		-- Query file inserted
    GROUP BY 
    	win_range, win_slide
)
SELECT 
	C.class, 
    E.num_queries, 
    count(D.*) AS num_speed_results, 
    round(max(diff_speed)::numeric,3) AS max_diff_speed, 
    round(avg(diff_speed)::numeric,3) AS avg_diff_speed, 
    round(stddev(diff_speed)::numeric,3) AS stdev_diff_speed
FROM axes_classification C, speed_differences D, queries_per_class E
WHERE 
	C.qid = D.qid
	AND C.class = E.win_slide
GROUP BY 
	C.class, E.num_queries
ORDER BY 
	C.class;
*/
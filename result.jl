files = [
"PID/40.out",
"Param/1.out",
"Param/2.out",
"Param/3.out",
"Param/4.out",
"Param/5.out",
"Param/6.out",
"Param/7.out",
"Param/8.out",
"Param/9.out",
"Param/10.out",
"Param/11.out",
"Param/12.out",
"Global/abelS100.out",
"Global/2_1_10S100.out",
"Global/2_1_7S100.out",
"Global/2_1_8S100.out",
"Global/5_2_5S100.out",
"Global/5_3_2S100.out",
"Global/8_4_1S100.out",
"Global/hydroS100.out",
"Global/immunS100.out",
"Global/st_fp7aS100.out",
"Global/st_fp7bS100.out",
"Global/st_fp7cS100.out",
"Global/st_fp7dS100.out",
"Global/st_fp7eS100.out",
"Global/st_fp8S100.out",
"Global/st_m1S100.out",
"Global/st_m2S100.out",
"Global/st_rv2S100.out",
"Global/st_rv3S100.out",
"Global/st_rv7S100.out",
"Global/st_rv8S100.out"];

NF = length(files)
Time = 1e5*ones(NF)
Nodes = 1e5*ones(Int, NF)
Results = 1e10*ones(Int, NF)
for i = 1:NF
   file = files[i]
   print(file)
   f = open(file) 
   nl = 1
   nl_node = -1
   while !eof(f)
     x = readline(f)
     if contains(x, "Solution time:  ")
     	number_as_string = split(x," ")[end-1]
	number = parse(Float64, number_as_string)
	Time[i] = number
     end
     if contains(x, "solved nodes:  ")
        number_as_string = split(x," ")[end]
        number = parse(Int, number_as_string)
        Nodes[i] = number
	nl_node = nl
     end
     if nl == (nl_node + 1)
        number_as_string = split(x)[end-2]
        number = parse(Float64, number_as_string)
        Results[i] = number
     end
     nl = nl + 1;
   end
   close(f)
end


println("Time   ")
for i = 1:NF
    println(Time[i])
end


println("Nodes   ")
for i = 1:NF
    println(Nodes[i])
end



Results_best = [
8.732906e+01
1.27628E+03
1.50022E+02
3.39626E+00
3.56219E+02
6.41359E+01
1.83498E+02
3.90203E+01
7.26087E+01
3.85143E+01
2.72312E+01
1.40485E+01
2.24499E+02
10069.4057
3067506.8886
-418330.0
552170.0
-354240.0
162.2201
18.3899
434673800.0000
0
-35655.5300
-63655.3000
-873107.3000
-11748.2000
-376380.0
1469100
-48781890.0000
-98942800.0000
-9161.0760
-3682.6350
-14842.9200
-13263.8200]


for i = 1:NF    
    if abs(Results_best[i] - Results[i])>=0.01 && abs(Results_best[i] - Results[i])/abs(Results_best[i])>=0.01
       println(files[i], "  best   ", Results_best[i], "    SNGO    ", Results[i])
    end
end

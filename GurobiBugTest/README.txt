When I test dominance method used in Danzig-Wolfe column geneartion, I found that in some instances, the number of iterations between using Gurobi and dominance for pricing problem are different.

For example, in t20_06_data.txt.

(1). in iteration 35, dominance method can't find result as good as gurobi.
gurobi:    route: 0,3,17,7,8,13,19; obj=2.9537
dominance: route: 0,3,7,8,13,19;    obj=2.9304

The reason is that a triangel inequality is broken: d(3,8)+d(8,7)-d(3,7)<0, d(7,8)+d(8,3)-d(7,3)<0.
After fixing it by decrease d(3,7) and d(7,3) by the diff (so that < becomes = in above equations), both gurobi and dominance find the same obj, 3.0565, the route is (0,3,9,7,5,13,19).

(2). after fixing triangel inequality problem, the number of iterations is still different. So I also tried forcing optimality tolerance to be the minimum, i.e.,1e-9.

(3). after (2), number of iterations is still differnt, then I noticed that in updated result, in iteration 43, gurobi can't find result as good as dominance method. I tried setting up the optimal route from dominance method as the initial solution in gurobi, then gurobi returns the better result the same as dominance. I think maybe this is a bug in gurobi. So I decide to export data for pricing problem in iteration 43 and solve it seperately to test if it is a real bug.

itr 43
gurobi:    route: (0,17,10,6,2,18,13,19); obj=0.0338
dominance: route: (0,17,3,8,10,15,19);    obj=0.042


SOLUTION

When I dig deeper, I found that the coefficients of x variables in pricing problem are slightly different. Then I found that in Iteration 31, Gurobi and dominace has two different routes with the same obj and distance.

itr 31
gurobi:    route: (0,17,15,11,9,5,13,19); obj=0.2298; dist=19.957
dominance: route: (0,17,15,11,9,13,19);   obj=0.2298; dist=19.957
xCoeff(9,5)=-0.261900
xCoeff(5,13)=-0.463200
xCoeffv(9,13)=-0.725100 = xCoeff(9,5) + xCoeff(5,13)
d(9,5)=2.619000
d(5,13)=4.632000
d(9,13)=7.251000 = d(9,5) + d(5,13)

It leads to differnt master models when using gurobi and dominanace method, then different solution in iteration 43.


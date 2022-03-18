function [n, no, lo, fo] = getatomicstates(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% =======
% Returns all electron states of the form (n, l, f) for a given Z
% Note: sum(fo) == Z
%
% Inputs
% ======
% Z             : atomic number
%
% Output
% ======
% no(:), lo(:)  : quantum numbers "n" and "l"
% fo(:)         : occupancy of the (n, l) state


switch Z
    
    case (1)
        n = 1;
        
        no = 1;
        lo = 0;
        fo = 1;
        
    case (2)
        n = 1;
        
        no = 1;
        lo = 0;
        fo = 2;
        
    case (3)
        n = 2;
        
        no = [ 1, 2 ];
        lo = [ 0, 0 ];
        fo = [ 2, 1 ];
        
    case (4)
        n = 2;
        
        no = [ 1, 2 ];
        lo = [ 0, 0 ];
        fo = [ 2, 2 ];
        
    case (5)
        n = 3;
        
        no = [ 1, 2, 2 ];
        lo = [ 0, 0, 1 ];
        fo = [ 2, 2, 1 ];
        
    case (6)
        n = 3;
        
        no = [ 1, 2, 2 ];
        lo = [ 0, 0, 1 ];
        fo = [ 2, 2, 2 ];
        
    case (7)
        n = 3;
        
        no = [ 1, 2, 2 ];
        lo = [ 0, 0, 1 ];
        fo = [ 2, 2, 3 ];
        
    case (8)
        n = 3;
        
        no = [ 1, 2, 2 ];
        lo = [ 0, 0, 1 ];
        fo = [ 2, 2, 4 ];
        
    case (9)
        n = 3;
        
        no = [ 1, 2, 2 ];
        lo = [ 0, 0, 1 ];
        fo = [ 2, 2, 5 ];
        
    case (10)
        n = 3;
        
        no = [ 1, 2, 2 ];
        lo = [ 0, 0, 1 ];
        fo = [ 2, 2, 6 ];
        
    case (11)
        n = 4;
        
        no = [ 1, 2, 2, 3 ];
        lo = [ 0, 0, 1, 0 ];
        fo = [ 2, 2, 6, 1 ];
        
    case (12)
        n = 4;
        
        no = [ 1, 2, 2, 3 ];
        lo = [ 0, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2 ];
        
    case (13)
        n = 5;
        
        no = [ 1, 2, 2, 3, 3 ];
        lo = [ 0, 0, 1, 0, 1 ];
        fo = [ 2, 2, 6, 2, 1 ];
        
    case (14)
        n = 5;
        
        no = [ 1, 2, 2, 3, 3 ];
        lo = [ 0, 0, 1, 0, 1 ];
        fo = [ 2, 2, 6, 2, 2 ];
        
    case (15)
        n = 5;
        
        no = [ 1, 2, 2, 3, 3 ];
        lo = [ 0, 0, 1, 0, 1 ];
        fo = [ 2, 2, 6, 2, 3 ];
        
    case (16)
        n = 5;
        
        no = [ 1, 2, 2, 3, 3 ];
        lo = [ 0, 0, 1, 0, 1 ];
        fo = [ 2, 2, 6, 2, 4 ];
        
    case (17)
        n = 5;
        
        no = [ 1, 2, 2, 3, 3 ];
        lo = [ 0, 0, 1, 0, 1 ];
        fo = [ 2, 2, 6, 2, 5 ];
        
    case (18)
        n = 5;
        
        no = [ 1, 2, 2, 3, 3 ];
        lo = [ 0, 0, 1, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6 ];
        
    case (19)
        n = 6;
        
        no = [ 1, 2, 2, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 1 ];
        
    case (20)
        n = 6;
        
        no = [ 1, 2, 2, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 2 ];
        
    case (21)
        n = 7;
        
        no = [ 1, 2, 2, 3, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 1, 2 ];
        
    case (22)
        n = 7;
        
        no = [ 1, 2, 2, 3, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 2, 2 ];
        
    case (23)
        n = 7;
        
        no = [ 1, 2, 2, 3, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 3, 2 ];
        
    case (24)
        n = 7;
        
        no = [ 1, 2, 2, 3, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 5, 1 ];
        
    case (25)
        n = 7;
        
        no = [ 1, 2, 2, 3, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 5, 2 ];
        
    case (26)
        n = 7;
        
        no = [ 1, 2, 2, 3, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 6, 2 ];
        
    case (27)
        n = 7;
        
        no = [ 1, 2, 2, 3, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 7, 2 ];
        
    case (28)
        n = 7;
        
        no = [ 1, 2, 2, 3, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 8, 2 ];
        
    case (29)
        n = 7;
        
        no = [ 1, 2, 2, 3, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 1 ];
        
    case (30)
        n = 7;
        
        no = [ 1, 2, 2, 3, 3, 3, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2 ];
        
    case (31)
        n = 8;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 1 ];
        
    case (32)
        n = 8;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 2 ];
        
    case (33)
        n = 8;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 3 ];
        
    case (34)
        n = 8;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 4 ];
        
    case (35)
        n = 8;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 5 ];
        
    case (36)
        n = 8;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6 ];
        
    case (37)
        n = 9;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 1 ];
        
    case (38)
        n = 9;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 2 ];
        
    case (39)
        n = 10;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 1, 2 ];
        
    case (40)
        n = 10;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 2, 2 ];
        
    case (41)
        n = 10;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 4, 1 ];
        
    case (42)
        n = 10;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 5, 1 ];
        
    case (43)
        n = 10;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 5, 2 ];
        
    case (44)
        n = 10;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 7, 1 ];
        
    case (45)
        n = 10;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 8, 1 ];
        
    case (46)
        n = 9;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10 ];
        
    case (47)
        n = 10;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 1 ];
        
    case (48)
        n = 10;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2 ];
        
    case (49)
        n = 11;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 1 ];
        
    case (50)
        n = 11;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 2 ];
        
    case (51)
        n = 11;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 3 ];
        
    case (52)
        n = 11;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 4 ];
        
    case (53)
        n = 11;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 5 ];
        
    case (54)
        n = 11;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6 ];
        
    case (55)
        n = 12;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6, 1 ];
        
    case (56)
        n = 12;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6, 2 ];
        
    case (57)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6, 1, 2 ];
        
    case (58)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 1, 2, 6, 1, 2 ];
        
    case (59)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 3, 2, 6, 2 ];
        
    case (60)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 4, 2, 6, 2 ];
        
    case (61)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 5, 2, 6, 2 ];
        
    case (62)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 6, 2, 6, 2 ];
        
    case (63)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 7, 2, 6, 2 ];
        
    case (64)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 7, 2, 6, 1, 2 ];
        
    case (65)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 9, 2, 6, 2 ];
        
    case (66)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 10, 2, 6, 2 ];
        
    case (67)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 11, 2, 6, 2 ];
        
    case (68)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 12, 2, 6, 2 ];
        
    case (69)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 13, 2, 6, 2 ];
        
    case (70)
        n = 13;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 2 ];
        
    case (71)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 1, 2 ];
        
    case (72)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 2, 2 ];
        
    case (73)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 3, 2 ];
        
    case (74)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 4, 2 ];
        
    case (75)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 5, 2 ];
        
    case (76)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 6, 2 ];
        
    case (77)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 7, 2 ];
        
    case (78)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 9, 1 ];
        
    case (79)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 1 ];
        
    case (80)
        n = 14;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2 ];
        
    case (81)
        n = 15;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 1 ];
        
    case (82)
        n = 15;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 2 ];
        
    case (83)
        n = 15;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 3 ];
        
    case (84)
        n = 15;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 4 ];
        
    case (85)
        n = 15;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 5 ];
        
    case (86)
        n = 15;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6 ];
        
    case (87)
        n = 16;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 1 ];
        
    case (88)
        n = 16;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 2 ];
        
    case (89)
        n = 17;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 1, 2 ];
        
    case (90)
        n = 17;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 2, 2 ];
        
    case (91)
        n = 18;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 2, 6, 1, 2 ];
        
    case (92)
        n = 18;
        
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7 ];
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0 ];
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 3, 2, 6, 1, 2 ];
        
    otherwise
        fprintf('Z = %-3d is an invalid Z-value\n', Z);
        error('Z not supported. Aborting . . .')
end

end
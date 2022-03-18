function setplotcommand(x,phi)

[m n] = size(phi);
s = [];
for i = 1:n
   if (i == 1)
     s0 = sprintf('x,phi(:,%d),',i);
   else
     s0 = sprintf(',x,phi(:,%d),',i);
   end
   if (i == 1)
     s1 = sprintf('''-b''');
   elseif (i == 2)
     s1 = sprintf('''-g''');
   elseif (i == 3)
     s1 = sprintf('''-r''');
   elseif (i == 4)
     s1 = sprintf('''--b''');
   elseif (i == 5)
     s1 = sprintf('''--g''');
   elseif (i == 6)
     s1 = sprintf('''--r''');
   elseif (i == 7)
     s1 = sprintf('''-ro''');
   elseif (i == 8)
     s1 = sprintf('''--ro''');
   elseif (i == 9)
     s1 = sprintf('''-go''');
   elseif (i == 10)
     s1 = sprintf('''--go''');
   elseif (i == 11)
     s1 = sprintf('''-rs''');
   elseif (i == 12)
     s1 = sprintf('''--rs''');
   elseif (i == 13)
     s1 = sprintf('''-gs''');
   elseif (i == 14)
     s1 = sprintf('''--gs''');
   else
     s1 = sprintf('''--''');
   end
   s = [s s0 s1];
end
eval(sprintf('plot(%s)',s));


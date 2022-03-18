function setlegend(n,type)

s = [];
for i = 1:n
   if (i < n)
      if (strcmp(type,'basis'))
        s0 = sprintf('''\\phi_{%d}'',',i);
      end
   else
      if (strcmp(type,'basis'))
        s0 = sprintf('''\\phi_{%d}''',i);
      end
   end
   s = [s s0];
end
eval(sprintf('legend(%s)',s));

return
end

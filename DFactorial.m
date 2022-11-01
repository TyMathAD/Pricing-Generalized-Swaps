function DFact = DFactorial(n)
  DFact = 1;
  if n > 1
    start = 1 + mod(n + 1,2); % caters for either parity of n
    DFact = prod(start:2:n);
  end
end
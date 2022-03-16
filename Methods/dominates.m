%Function to check if x dominates y
function answer = dominates(Solution1,Solution2)
answer = all(Solution1.y <= Solution2.y) && any(Solution1.y<Solution2.y);
end
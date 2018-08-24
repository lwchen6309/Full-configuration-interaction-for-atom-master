function Combination = Ball2Box(BallNum,BoxNum)
if BoxNum > 0
    SizeOfCombination = factorial(BallNum + BoxNum - 1) / ...
    factorial(BallNum) / factorial(BoxNum - 1);
    Combination = zeros(1,BoxNum);
    while BallNum > 0
        NewCombination = [];
        for i = Combination.'
            % Add one ball to List and push to NewCombination
            NewCombination = [NewCombination,AddOneBall(i)];
        end
        Combination = NewCombination.';
        BallNum = BallNum - 1;
    end
    Combination = flipud(unique(Combination,'rows'));
end
assert(SizeOfCombination == size(Combination,1));
end

function NewCombination = AddOneBall(Combination)
    SizeOfBox = size(Combination,1);
    NewCombination = repmat(Combination,1,SizeOfBox);
    NewCombination = NewCombination + eye(SizeOfBox);
end
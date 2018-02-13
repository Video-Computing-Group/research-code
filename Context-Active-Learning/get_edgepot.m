function newEdgePot = get_edgepot(edgepot, Link, Index, Label)

newEdgePot = zeros(size(edgepot));

for i = 1:size(Link,1)
    t1 = Link(i,1)==Index;
    t2 = Link(i,2)==Index;
    if sum(t1) && sum(t2)
        newEdgePot(Label(t2),Label(t1)) = newEdgePot(Label(t2),Label(t1)) + 1;
    end
end

newEdgePot = newEdgePot + newEdgePot';
newEdgePot = newEdgePot - 0.5*eye(size(newEdgePot)).*newEdgePot;  % Coz. diag elements has been added twice

newEdgePot = newEdgePot + edgepot;
function draw_solution(nodes,edges,u)
figure
hold on

trimesh(edges',nodes(1,:),nodes(2,:),u)

view(-35,20)
xlabel('x_1')
ylabel('x_2')
zlabel('u(x_1,x_2)')
hold off
end
function [archive_pos, archive_val] = update_archive(archive_positions, archive_values)
    % Non-dominated sorting to update the archive
    non_dominated = true(size(archive_values, 1), 1);
    for i = 1:size(archive_values, 1)
        for j = 1:size(archive_values, 1)
            if i ~= j && dominates(archive_values(j, :), archive_values(i, :))
                non_dominated(i) = false;
                break;
            end
        end
    end
    archive_pos = archive_positions(non_dominated, :);
    archive_val = archive_values(non_dominated, :);
end
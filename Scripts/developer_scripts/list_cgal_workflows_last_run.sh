#!/bin/bash
echo "| repo |  workflow  | branch | event | runs on  |  status of last run  | state | annotation |  date  | date since last runs |  file  |"
echo "| :--: | :--------: | :----: | :---: | :------: | :------------------: | :---: | :--------: | :----: | :------------------: | :----: |"
actualdate=$EPOCHSECONDS
for repo in $(gh api orgs/CGAL/repos --jq '.[].full_name' | grep -v dev )
do
    if [ "$repo" != "CGAL/CNRS" ] && [ "$repo" != "CGAL/GeometryFactory" ]
    then
        default_branch=$(gh api repos/$repo --jq '.default_branch')
        workflows=$(gh api repos/$repo/actions/workflows)
        workflows_count=$(jq '.total_count' <<< "$workflows")
        for ((i=0;i<workflows_count;i++))
        do
            workflow_id=$(jq '.workflows['$i'].id' <<< "$workflows")
            workflows_state=$(jq '.workflows['$i'].state' <<< "$workflows")
            workflow_runs=$(gh api repos/$repo/actions/workflows/$workflow_id/runs)
            workflows_name=$(jq -e '.workflow_runs[0].name' <<< "$workflow_runs")
            if [ $? -eq 0 ]
            then
                workflows_status=$(jq -r '.workflow_runs[0].status' <<< "$workflow_runs")
                workflows_conclusion=$(jq -r '.workflow_runs[0].conclusion' <<< "$workflow_runs")
                workflows_start=$(jq -r '.workflow_runs[0].run_started_at' <<< "$workflow_runs")
                workflows_date=$( date --date="$workflows_start" +%s )
                workflows_on=$(jq -r '.workflow_runs[0].event' <<< "$workflow_runs")
                workflows_path=$(jq -r '.workflow_runs[0].path' <<< "$workflow_runs")
                workflows_branch=$(jq -r '.workflow_runs[0].head_branch' <<< "$workflow_runs")
                workflows_checksuite_id=$(jq -r '.workflow_runs[0].check_suite_id' <<< "$workflow_runs")
                workflows_check_runs=$(gh api repos/$repo/check-suites/$workflows_checksuite_id/check-runs)
                workflows_check_runs_id=$(jq -r '.check_runs[0].id' <<< "$workflows_check_runs")
                workflows_check_runs_annotation=$(gh api repos/$repo/check-runs/$workflows_check_runs_id/annotations)
                worfklows_annotation_level=$(jq -r '.[].annotation_level' <<< "$workflows_check_runs_annotation")
                if [ "$worfklows_annotation_level" == "" ]
                then
                    worfklows_annotation_level+="-"
                fi
                workflows_event=""
                for trigger in $(curl --silent https://raw.githubusercontent.com/$repo/$default_branch/$workflows_path | yq '.on' | grep -v '  .*'); do
                    if [ "${trigger}" != "-" ]
                    then
                        workflows_event+=${trigger}
                        workflows_event+="<br/>"
                    fi
                done
                echo "| $repo | $workflows_name | $workflows_branch | $workflows_event | $workflows_on | $workflows_status - $workflows_conclusion | $workflows_state | ***$worfklows_annotation_level*** | $workflows_start | $(((actualdate - workflows_date) / 86400 )) days | $repo/$workflows_path"
            fi
        done
    fi
done

echo ""
echo "| repo | dependabot |"
echo "| :--: | :--------: |"
for repo in $(gh api orgs/CGAL/repos --jq '.[].full_name' | grep -v dev )
do
    if [ "$repo" != "CGAL/CNRS" ] && [ "$repo" != "CGAL/GeometryFactory" ]
    then
        default_branch=$(gh api repos/$repo --jq '.default_branch')
        workflows=$(gh api repos/$repo/actions/workflows)
        workflows_count=$(jq '.total_count' <<< "$workflows")
        if [ $workflows_count != 0 ]
        then
            dependabot=$(curl --silent https://raw.githubusercontent.com/$repo/$default_branch/.github/dependabot.yml)
            dependabotexist=""
            if [ "$dependabot" != "404: Not Found" ]
            then
                dependabotexist="yes"
            else
                dependabotexist="no"
            fi
            echo "| $repo | $dependabotexist |"
        fi
    fi
done

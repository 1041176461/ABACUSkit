'''
Date: 2021-02-14 09:07:04
LastEditors: jiyuyang
LastEditTime: 2021-03-05 17:40:10
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import abc
from collections import defaultdict


class JobResource(abc.ABC):
    """Data structure to store job resources.

    Typical attributes are:
    * ``num_machines``
    * ``num_mpiprocs_per_machine``

    or (e.g. for SGE)
    * ``tot_num_mpiprocs``
    * ``parallel_env``

    The constructor should take care of checking the values.
    """

    _default_fields = tuple()

    @abc.abstractmethod
    def validate_resources(cls, **kwargs):
        """Validate the resources against the job resource class of this scheduler.

        :param kwargs: dictionary of values to define the job resources
        :raises ValueError: if the resources are invalid or incomplete
        :return: optional tuple of parsed resource settings
        """
        resources = defaultdict(lambda: None)

        return resources


class NodeNumberJobResource(JobResource):
    """`JobResource` for schedulers that support the specification of a number of nodes and cpus per node."""

    _default_fields = (
        'num_machines',
        'num_mpiprocs_per_machine',
        'num_cores_per_machine',
        'num_cores_per_mpiproc',
    )

    @classmethod
    def validate_resources(cls, **kwargs):
        """Validate the resources against the job resource class of this scheduler.

        :param kwargs: dictionary of values to define the job resources
        :return: defaultdict with the parsed parameters populated
        :raises ValueError: if the resources are invalid or incomplete
        """
        resources = defaultdict(lambda: None)

        def is_greater_equal_one(parameter):
            value = resources[parameter]
            if value is not None and value < 1:
                raise ValueError(
                    f'`{parameter}` must be greater than or equal to one.')

        for parameter in list(cls._default_fields) + ['tot_num_mpiprocs']:
            value = kwargs.pop(parameter, None)
            if value == None:
                resources[parameter] = value
            else:
                try:
                    resources[parameter] = int(value)
                except ValueError:
                    raise ValueError(
                        f'`{parameter}` must be an integer when specified')

        if kwargs:
            raise ValueError(
                f"these parameters were not recognized: {','.join(list(kwargs.keys()))}")

        if [resources['num_machines'], resources['num_mpiprocs_per_machine'], resources['tot_num_mpiprocs']].count(None) > 1:
            raise ValueError(
                'At least two among `num_machines`, `num_mpiprocs_per_machine` or `tot_num_mpiprocs` must be specified'
            )

        for parameter in ['num_machines', 'num_mpiprocs_per_machine']:
            is_greater_equal_one(parameter)

        if resources['num_machines'] is None:
            resources['num_machines'] = resources['tot_num_mpiprocs'] // resources['num_mpiprocs_per_machine']
        elif resources['num_mpiprocs_per_machine'] is None:
            resources['num_mpiprocs_per_machine'] = resources['tot_num_mpiprocs'] // resources['num_machines']
        elif resources['tot_num_mpiprocs'] is None:
            resources['tot_num_mpiprocs'] = resources['num_mpiprocs_per_machine'] * \
                resources['num_machines']

        if resources['tot_num_mpiprocs'] != resources['num_mpiprocs_per_machine'] * resources['num_machines']:
            raise ValueError(
                '`tot_num_mpiprocs` is not equal to `num_mpiprocs_per_machine * num_machines`.')

        is_greater_equal_one('num_mpiprocs_per_machine')
        is_greater_equal_one('num_machines')

        return resources

    def __init__(self, **kwargs):
        """Initialize the job resources from the passed arguments.

        :raises ValueError: if the resources are invalid or incomplete
        """
        self.resources = self.validate_resources(**kwargs)


class ParEnvJobResource(JobResource):
    """`JobResource` for schedulers that support the specification of a parallel environment and number of MPI procs."""

    _default_fields = (
        'parallel_env',
        'tot_num_mpiprocs',
    )

    @classmethod
    def validate_resources(cls, **kwargs):
        """Validate the resources against the job resource class of this scheduler.

        :param kwargs: dictionary of values to define the job resources
        :return: defaultdict with the parsed parameters populated
        :raises ValueError: if the resources are invalid or incomplete
        """
        resources = defaultdict(lambda: None)

        try:
            resources["parallel_env"] = kwargs.pop('parallel_env')
        except KeyError:
            raise ValueError(
                '`parallel_env` must be specified and must be a string')
        else:
            if not isinstance(resources["parallel_env"], str):
                raise ValueError(
                    '`parallel_env` must be specified and must be a string')

        try:
            resources["tot_num_mpiprocs"] = int(kwargs.pop('tot_num_mpiprocs'))
        except (KeyError, TypeError, ValueError):
            raise ValueError(
                '`tot_num_mpiprocs` must be specified and must be an integer')

        if resources["tot_num_mpiprocs"] < 1:
            raise ValueError(
                '`tot_num_mpiprocs` must be greater than or equal to one.')

        if kwargs:
            raise ValueError(
                f"these parameters were not recognized: {', '.join(list(kwargs.keys()))}")

        return resources

    def __init__(self, **kwargs):
        """Initialize the job resources from the passed arguments.

        :raises ValueError: if the resources are invalid or incomplete
        """
        self.resources = self.validate_resources(**kwargs)
